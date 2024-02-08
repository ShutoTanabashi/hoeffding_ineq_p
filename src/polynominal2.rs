//! 有薗ら提案(2022)
//! Markov近似の際に多項式を分母に用いるHoeffding型確立不等式．n=2に限定した版．

// テスト用
use super::hoeffding::HoeffdingIneq;

use super::{IneqError, UpperBound, CompareHoeffding};

extern crate libm;
use libm::{log};
extern crate rayon;
use rayon::prelude::*;
extern crate simple_excel_writer;
use simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};

/// 方程式の解 $Z_1, Z_2$ を示す型  
/// 次の3種類のうち，どの分類であるかを示す．
/// 
/// - 重根 `Multiple`
/// - 複号同順根 `Real`
/// - エッジ根（値域の端） `Edge`
/// 
/// いずれにおいても解の値を `f64` 型で保持する．
#[derive(Debug, Clone, Copy)]
pub enum Root {
    Multiple(f64),
    Real(f64, f64),
    Edge(f64, f64),
}


impl Root {
    /// 解の値をタプルで返す
    pub fn value(&self) -> (f64, f64) {
        match self {
            Root::Multiple(z) => (z.clone(), z.clone()),
            Root::Real(z1, z2) => Root::ord_z(z1, z2),
            Root::Edge(z1, z2) => Root::ord_z(z1, z2),
        }
    }

    // 2個の引数について，値の大小を並び変えたタプルを返す
    // self.value用の内部関数．
    fn ord_z(z1: &f64, z2: &f64) -> (f64, f64) {
        if z1 < z2 {
            (z1.clone(), z2.clone())
        } else {
            (z2.clone(), z1.clone())
        }
    }

    /// 解の分類を日本語のStringを出力
    pub fn kind_root(&self) -> String {
        match self {
            Root::Multiple(_) => "重根".to_string(),
            Root::Real(_,_) => "複号異順解".to_string(),
            Root::Edge(_,_) => "エッジ解".to_string(),
        }
    }

}


// PolynominalIneqの内部計算用
// 計算結果とそれに用いた内部変数をセットで格納する
#[allow(dead_code)]
struct VerbosePolynominalIneq {
    // 計算対象のz
    z: f64,
    // 上側確率の上界
    pr: f64,
    // 関数 $\eta(Z)$ の停留点
    root: Root,
    // 関数 $\eta(Z)$ の係数
    coef_eta: (f64,f64,f64),
    // 定数$\theta$
    theta: f64,
    // 計算に用いた$\prod \eta(Z_i)$ の値
    pi_eta: f64,
}

impl VerbosePolynominalIneq {
    // 関数 $\eta$の2次の項の係数 $alpha$
    #[allow(dead_code)]
    pub fn alpha(&self) -> f64 {
        self.coef_eta.0
    }

    // 関数 $\eta$の1次の項の係数 $beta$
    #[allow(dead_code)]
    pub fn beta(&self) -> f64 {
        self.coef_eta.1
    }

    // 関数 $\eta$の0次の項の係数 $gamma$
    #[allow(dead_code)]
    pub fn gamma(&self) -> f64 {
        self.coef_eta.2
    }

    // $\eta(Z)$ の停留点の判定式に用いる値 $2 \alpha Z_i^{\dag} + \beta$ の値を計算
    #[allow(dead_code)]
    fn calc_descriminat(&self, z: &f64) -> f64 {
        let alpha = self.alpha();
        let beta = self.beta();
        2.0 * alpha * z + beta
    }

    // $\eta(Z)$ の停留点の判定式に用いる値 $2 \alpha Z_i^{\dag} + \beta$ の値をタプルで返す
    #[allow(dead_code)]
    pub fn discriminant(&self) -> (f64, f64) {
        let (z1, z2) = self.root.value();
        (self.calc_descriminat(&z1), self.calc_descriminat(&z2))
    }
}

/// 有薗らによる多項式を分母に用いた確立不等式  
/// 現時点では $n=2$ の場合にのみ対応．
#[derive(Debug, Clone)]
pub struct PolynominalIneq {
    // 確率変数の個数
    n: usize,
    // 確率変数の分散
    sigma2: f64,
    // 確率変数の最大値
    delta: f64,
}

impl PolynominalIneq {
    /// 各種引数から有薗らによる多項式を分母に用いた確立不等式を初期化  
    /// $n=2$版．
    /// 
    /// # 引数
    /// * `sigma2` - 確率変数の分散の平均
    /// * `delta` - 確率変数の定義域の最大値
    /// 
    /// # 使用例
    /// ```
    /// # use hoeffding_ineq::polynominal2::PolynominalIneq;
    /// let polynominal = PolynominalIneq::new(4.0, 5.0).unwrap();
    /// ```
    pub fn new(sigma2: f64, delta: f64) -> Result<Self, IneqError> {
        if sigma2 <= 0.0 {
            Err(IneqError{
                message: "sigma2 must be greater than zero.".to_string()
            })
        } else {
            Ok(PolynominalIneq{n: 2, sigma2, delta})
        }
    }


    // Hoeffdingの確立不等式にてネイピア数 $e$ の指数にある定数 $\theta$ を計算
    // 
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    //
    // # 定義式(TeX)
    // > \theta = \frac{n \delta}{\delta^2 + \sigma^2} 
    // > \times 
    // > \ln{\frac{1 + \frac{z \delta}{\sigma^2}}{1-\frac{z}{\delta}}}
    fn theta(&self, z: &f64) -> f64 {
        (self.n as f64) * self.delta / (self.delta.powi(2) + self.sigma2) 
        * log( (1.0 + z * self.delta / self.sigma2) / (1.0 - z / self.delta))
    } 
    

    // 確立不等式の分母となる多項式関数 $\eta(Z)$ 内の3個の係数 $\alpha, \beta, \gamma$ を計算
    //
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    //
    // # 返り値
    // * `(alpha, beta, gamma)` - 下記の $\eta(Z)$ の係数
    //
    // # 定義式(TeX)
    // > \eta(Z) = \alpha Z^2 + \beta Z + \gamma
    //
    // > \alpha = \frac{\delta^2}{(\delta^2 + \sigma^2)^2}
    // >          (e^{\delta \frac{\theta}{n}} 
    // >                - e^{- \frac{\sigma^2}{\delta} \frac{\theta}{n} })
    // >          -
    // >          \frac{\delta}{\delta^2 + \sigma^2} \frac{\theta}{n}
    // >                e^{- \frac{\sigma^2}{\delta} \frac{\theta}{n} }
    //
    // > \beta = \frac{2 \delta \sigma^2}{(\delta^2 + \sigma^2)^2}
    // >            (e^{\delta \frac{\theta}{n}} 
    // >                - e^{- \frac{\sigma^2}{\delta} \frac{\theta}{n} })
    // >         +
    // >         \frac{\delta^2 - \sigma^2}{\delta^2 + \sigma^2} 
    // >            \frac{\theta}{n}
    // >            e^{- \frac{\sigma^2}{\delta} \frac{\theta}{n} }
    //
    // > \gamma = \frac{(\sigma^2)^2}{(\delta^2 + \sigma^2)^2}
    // >         e^{\delta \frac{\theta}{n}} 
    // >         +
    // >         (\frac{\delta^2 + 2 \sigma^2}{(\delta^2 + \sigma^2)^2}
    // >            + \frac{\delta \sigma^2}{\delta^2 + \sigma^2}
    // >                \frac{\theta}{n})
    // >            e^{- \frac{\sigma^2}{\delta} \frac{\theta}{n} }
    fn cal_coef_eta(&self, z: &f64) -> (f64, f64, f64) {
        let theta = self.theta(z);
        let d2 = self.delta.powi(2);
        let d2ps2 = d2 + self.sigma2;
        let sqd2ps2 = d2ps2.powi(2);
        let tdivn = theta / (self.n as f64);
        let edtn = (self.delta * tdivn).exp();
        let estdn = (- tdivn * self.sigma2 / self.delta).exp();

        let alpha = d2 / sqd2ps2 * (edtn - estdn)
                    -
                    self.delta / d2ps2 * tdivn * estdn;
        let beta = 2.0 * self.delta * self.sigma2 / sqd2ps2  * (edtn - estdn)
                   +
                   (d2 - self.sigma2) / d2ps2 * tdivn * estdn;
        let gamma = self.sigma2.powi(2) / sqd2ps2 * edtn
                    +
                    (d2 * (d2 + 2.0 * self.sigma2) / sqd2ps2
                     +
                     self.delta * self.sigma2 / d2ps2 * tdivn)
                    *
                    estdn;
        
        (alpha, beta, gamma)
    }


    // 関数 $\eta(Z)$ の停留点を与える解 $Z_1, Z_2$ の計算
    //
    // 内部計算用の関数であり，引数に $\eta(Z)$ の係数 $\alpha, \beta, \gamma$ を取る．
    //
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    //
    // # 定義式(TeX)
    // ## 実数解
    // 条件: 解を与える式
    // 
    // > (2 \aplha Z_i + \beta) = (2 \alpha \bar{Z} + \beta) 
    // >                          \pm
    // >                          \sqrt{(2 \alpha \bar{Z} + \beta)^2 - (4 \alpha \gamma - \beta^2)}
    // 
    // が実数解を持ち，大きい方の解が確率変数の値域 $\delta$ を超えない場合．
    //
    // ## エッジ解
    // 条件: 解を与える式が実数解を持ち，大きい方の解が確率変数の値域 $\delta$ を超える場合．
    // この場合は次の根が与えられる．
    //
    // > (Z_1, Z_2) = (\delta, 2 \bar{Z} - \delta)
    //
    // ## 重根
    // 条件 解を与える式が虚数解を持つ場合．
    // この場合は次の式から根が与えられる．
    //
    // > (2 \aplha Z_i + \beta) = (2 \alpha \bar{Z} + \beta) 
    fn calc_root_stationary(&self, z: &f64, coef_eta: &(f64, f64, f64)) -> Root {
        let (alpha, beta, gamma) = coef_eta;
        let root_mul_cnd = 2.0 * alpha * z + beta; 
        let discri = root_mul_cnd.powi(2) - (4.0 * alpha * gamma - beta.powi(2));

        if discri < 0.0 {
            Root::Multiple( z.clone() )
        } else {
            let threshold = root_mul_cnd + discri.sqrt();
            if threshold > (2.0 * alpha * self.delta + beta) {
                Root::Edge(self.delta, 2.0 * z - self.delta)
            } else {
                let root_real_max = z + discri.sqrt() / (2.0 * alpha);
                let root_real_min = z - discri.sqrt() / (2.0 * alpha);
                Root::Real(root_real_min, root_real_max)
            }
        }
    }

    /// 与えられた変数 $z$ に対して，確立不等式の分母となる関数 $\eta(Z)$ の停留点を与える解
    /// 
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn root_stationary(&self, z: &f64) -> Root {
        let coef_eta = self.cal_coef_eta(z);
        self.calc_root_stationary(z, &coef_eta)
    }


    // 確立不等式の分母となる，停留点での $\eta(Z_i)$ の値の積
    // 内部計算用のため，すべての変数を引数に取ります．
    // 利用する際は適宜ラッパー関数を利用してください．
    // 
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値 
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    // * `root` - 停留点での $Z_1, Z_2$
    fn calc_pi_eta(&self, z: &f64, coef_eta: &(f64, f64, f64), root: &Root) -> f64 {
        let (alpha, beta, gamma) = coef_eta;
        match root {
            Root::Multiple(_) => {
                (
                    (2.0 * alpha * z + beta).powi(2)
                    +
                    4.0 * alpha * gamma - beta.powi(2)
                ).powi(2)
                / ((4.0 * alpha).powi(2))
            },
            Root::Real(_,_) => {
                (2.0 * alpha * z + beta).powi(2)
                *
                (4.0 * alpha * gamma - beta.powi(2))
                / ((2.0 * alpha).powi(2))
            },
            Root::Edge(_,_) => {
                (
                    (2.0 * alpha * self.delta + beta).powi(2)
                    +
                    4.0 * alpha * gamma - beta.powi(2)
                ) * (
                    (2.0 * alpha * (2.0 * z - self.delta) + beta).powi(2)
                    +
                    (4.0 * alpha * gamma - beta.powi(2))
                )
                / ((4.0 * alpha).powi(2))
            }
        }
    }

    // 上側確率の上界を計算
    // 内部計算用のため，すべての変数を引数に取ります．
    // 利用する際は適宜ラッパー関数を利用してください．
    // 
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値 
    fn calculate(&self, z: &f64) -> (f64, Root) {
        let coef_eta = self.cal_coef_eta(z);
        let root = self.calc_root_stationary(z, &coef_eta);
        let pi_eta = self.calc_pi_eta(z, &coef_eta, &root);
        let theta = self.theta(z);
        let d2ps2 = self.delta.powi(2) + self.sigma2;

        let prob = ((
            self.delta.powi(2) / d2ps2 * ((- theta / (self.n as f64)* self.sigma2 / self.delta).exp())
            +
            self.sigma2 / d2ps2 * (( theta / (self.n as f64) * self.delta).exp())
            ).powi(self.n as i32)) / pi_eta;

        (prob, root)
    }


    /// 上側確率の上界と計算に用いた $\eta(Z_i)$ の根を取得
    /// 
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn pr_with_root(&self, z: &f64) -> Result<(f64, Root), IneqError> {
        if *z <= self.min_z() {
            Err(IneqError {
                message: format!("`z` must be greater than {}", self.min_z())
            })
        } else if *z >= self.max_z() {
            Err(IneqError {
                message: format!("`z` must be less than {}", self.max_z())
            })
        } else {
            Ok(self.calculate(z))
        }
    }


    /// 複数個の変数に対して上側確率の上界と計算に用いた $\eta(Z_i)$ の根を計算
    /// 
    /// # 引数
    /// * `zs` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値の組
    pub fn map_pr_root(&self, zs: &[f64]) -> Result<Vec<(f64, Root)>, IneqError> {
        zs.par_iter()
          .map(|z| self.pr_with_root(z))
          .collect::<Result<Vec<(f64, Root)>, IneqError>>()
    }


    /// 上側確率の上界を計算し，(z, Pr(z), Root)のタプルを返す
    /// 
    /// # 引数  
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn pr_with_root_tuple(&self, z: &f64) -> Result<(f64, f64, Root), IneqError> {
        let (pr, root) = self.pr_with_root(z)?;
        Ok((z.clone(), pr, root))
    }


    /// 引数の有効な範囲内における上側確率の上界をまとめて計算
    /// 返り値のタプルは(z, Pr(z), Root)を意味する
    pub fn overview_with_root(&self) -> Result<Vec<(f64, f64, Root)>, IneqError> {
        let vec_z = self.variable_for_overview();
        vec_z.par_iter()
             .map(|z| self.pr_with_root_tuple(z))
             .collect::<Result<Vec<(f64, f64, Root)>, IneqError>>()
    }


    // 引数z に対してVerbosePolynominalIneqを計算
    //
    // # 引数  
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    fn pr_with_verbose(&self, z: &f64) -> Result<VerbosePolynominalIneq, IneqError> {
        let (pr, root) = self.pr_with_root(z)?;
        let coef_eta = self.cal_coef_eta(z);
        let pi_eta = self.calc_pi_eta(z, &coef_eta, &root);
        let theta = self.theta(z);
        Ok(VerbosePolynominalIneq{z: *z, pr, root, coef_eta, pi_eta, theta})
    }


    // 引数の有効な範囲内における上側確率の上界を，詳細情報とともにまとめて計算  
    fn overview_with_verbose(&self) -> Result<Vec<VerbosePolynominalIneq>, IneqError> {
        let vec_z = self.variable_for_overview();
        vec_z.par_iter()
             . map(|z| self.pr_with_verbose(z))
             .collect::<Result<Vec<VerbosePolynominalIneq>, IneqError>>()
    }

    /// 引数の有効な範囲内における上側確率の上界と計算に用いた $\eta(Z_i)$ の根をExcelファイルで保存
    /// 
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    pub fn overview_with_root_to_excel(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let overview = self.overview_with_verbose()?;
        let params = self.param_to_tuple();


        // Excelへの書き込み

        let mut wb = xlsx::Workbook::create(xlsx_path);
        let mut sheet_1 = wb.create_sheet("Calculated Value");
        
        wb.write_sheet(&mut sheet_1, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "提案手法","2αZ1+β", "2αZ2+β", "Z1", "Z2", "根の種類"])?;
            for ver in overview.iter() {
                let (z1, z2) = ver.root.value();
                let (cal_z1, cal_z2) = ver.discriminant();
                sheet_writer.append_row(xlsx::row![
                        ver.z.to_cell_value(), 
                        ver.pr.to_cell_value(), 
                        cal_z1.to_cell_value(),
                        cal_z2.to_cell_value(),
                        z1.to_cell_value(),
                        z2.to_cell_value(),
                        ver.root.kind_root()
                ])?;
            }
            Ok(())
        })?;

        let mut sheet_2 = wb.create_sheet("Parameters");
        
        wb.write_sheet(&mut sheet_2, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["Parameter", "Value"])?;
            for (p, v) in params.iter() {
                sheet_writer.append_row(xlsx::row![p.to_cell_value(), v.to_cell_value()])?;
            }
            Ok(())
        })?;
        

        // 計算に用いた内部変数
        let mut sheet_3 = wb.create_sheet("Internal value");
        
        wb.write_sheet(&mut sheet_3, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Z1", "Z2", "2αZ1+β", "2αZ2+β", "種類", "theta", "alpha", "beta", "gamma", "Pi eta(Z_i)"])?;
            for ver in overview.iter() {
                let (z1, z2) = ver.root.value();
                let (cal_z1, cal_z2) = ver.discriminant();
                sheet_writer.append_row(xlsx::row![
                    ver.z.to_cell_value(), 
                    ver.pr.to_cell_value(),
                    z1.to_cell_value(),
                    z2.to_cell_value(),
                    cal_z1.to_cell_value(),
                    cal_z2.to_cell_value(),
                    ver.root.kind_root(),
                    ver.theta.to_cell_value(),
                    ver.alpha().to_cell_value(), 
                    ver.beta().to_cell_value(), 
                    ver.gamma().to_cell_value(),
                    ver.pi_eta.to_cell_value()
                    ])?;
            }
            Ok(())
        })?;


        Ok(())
        
    }


    /// 提案手法とHoeffdingの確率不等式を比較し，Excelで保存
    ///
    /// # 引数
    /// `xlsx_path` - 保存先のExcelファイルパス
    pub fn compare_hoeffding_to_excel(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let overview = self.overview_with_verbose()?;
        let params = self.param_to_tuple();


        // Excelへの書き込み

        let mut wb = xlsx::Workbook::create(xlsx_path);
        let mut sheet_1 = wb.create_sheet("Calculated Value");
        
        let hoeff = self.same_condition_hoeffding()?;
        wb.write_sheet(&mut sheet_1, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Hoeffding","2αZ1+β", "2αZ2+β", "根の種類"])?;
            for ver in overview.iter() {
                let (cal_z1, cal_z2) = ver.discriminant();
                let pr_hoeff = hoeff.pr(&ver.z);
                sheet_writer.append_row(xlsx::row![
                        ver.z.to_cell_value(), 
                        ver.pr.to_cell_value(), 
                        pr_hoeff.to_cell_value(),
                        cal_z1.to_cell_value(),
                        cal_z2.to_cell_value(),
                        ver.root.kind_root()
                ])?;
            }
            Ok(())
        })?;

        let mut sheet_2 = wb.create_sheet("Parameters");
        
        wb.write_sheet(&mut sheet_2, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["Parameter", "Value"])?;
            for (p, v) in params.iter() {
                sheet_writer.append_row(xlsx::row![p.to_cell_value(), v.to_cell_value()])?;
            }
            Ok(())
        })?;
        

        // 計算に用いた内部変数
        let mut sheet_3 = wb.create_sheet("Internal value");
        
        wb.write_sheet(&mut sheet_3, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Z1", "Z2", "2αZ1+β", "2αZ2+β", "種類", "theta", "alpha", "beta", "gamma", "Pi eta(Z_i)"])?;
            for ver in overview.iter() {
                let (z1, z2) = ver.root.value();
                let (cal_z1, cal_z2) = ver.discriminant();
                sheet_writer.append_row(xlsx::row![
                    ver.z.to_cell_value(), 
                    ver.pr.to_cell_value(),
                    z1.to_cell_value(),
                    z2.to_cell_value(),
                    cal_z1.to_cell_value(),
                    cal_z2.to_cell_value(),
                    ver.root.kind_root(),
                    ver.theta.to_cell_value(),
                    ver.alpha().to_cell_value(), 
                    ver.beta().to_cell_value(), 
                    ver.gamma().to_cell_value(),
                    ver.pi_eta.to_cell_value()
                    ])?;
            }
            Ok(())
        })?;

        Ok(())
    }
}

impl UpperBound for PolynominalIneq {
    fn pr(&self, z: &f64) -> f64 {
        if *z <= self.min_z() {
            1.0
        } else if *z >= self.max_z() {
            0.0
        } else {
            let (prob, _) = self.calculate(z);
            prob
        }
    }

    fn max_z(&self) -> f64 {
        self.delta
    }

    fn min_z(&self) -> f64 {
        0.0
    }

    fn param_to_tuple(&self) -> Vec<(String, f64)> {
        vec![
            ("n".to_string(), self.n as f64),
            ("sigma^2".to_string(), self.sigma2),
            ("delta".to_string(), self.delta)
        ]
    }
}


impl CompareHoeffding for PolynominalIneq {
    fn same_condition_hoeffding(&self) -> Result<HoeffdingIneq, IneqError> {
        let n = self.n;
        let sigma2 = self.sigma2;
        let delta = self.delta;
        HoeffdingIneq::new(n, sigma2, delta)
    }
}
