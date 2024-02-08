//! 有薗ら提案(2022)
//! Markov近似の際に多項式を分母に用いるHoeffding型確立不等式

// テスト用
use super::hoeffding::HoeffdingIneq;

use super::{IneqError, UpperBound, CompareHoeffding};


extern crate libm;
use libm::log;
extern crate rayon;
use rayon::prelude::*;
extern crate simple_excel_writer;
use simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};

/// 方程式の解$ (2 \aplha Z_i^{\ast} + \beta) $を示す型  
/// 次の3種類のうち，どの分類であるかを示す．
/// 
/// - 重根 `Multiple`
/// - 複号根 `Real`
/// - エッジ根（値域の端） `Edge`
/// 
/// いずれにおいても解の個数を`usize`型，解の値を `f64` 型で保持する．
#[derive(Debug, Clone)]
pub enum Root {
    Multiple(usize, f64),
    Real((usize, f64), (usize, f64)),
    Edge(usize, f64, Box<Root>),
}


impl Root {
    /// 解の値を[`Vec`]で返す
    pub fn value(&self) -> Vec<f64> {
        match self {
            Root::Multiple(n,z) => {
                    Root::vec_z(n, z)
                },
            Root::Real(a1, a2) => {
                    let ((n1, z1), (n2, z2)) =  Root::ord_z(a1, a2);
                    let mut roots_z1 = Root::vec_z(&n1, &z1);
                    let mut roots_z2 = Root::vec_z(&n2, &z2);
                    roots_z1.append(&mut roots_z2);
                    roots_z1
                },
            Root::Edge(n, z, rest_roots) => {
                    let mut roots_z = Root::vec_z(n, z);
                    let mut roots_rest = rest_roots.value();
                    roots_z.append(&mut roots_rest);
                    roots_z
            },
        }
    }

    fn vec_z(n: &usize, z: &f64) -> Vec<f64> {
        (0..*n).map(|_| *z)
               .collect()
    }

    // 2個の引数について，値の大小を並び変えたタプルを返す
    // self.value用の内部関数．
    fn ord_z(a1: &(usize, f64), a2: &(usize, f64)) -> ((usize, f64), (usize, f64)) {
        if a1.1 > a2.1 {
            (a1.clone(), a2.clone())
        } else {
            (a2.clone(), a1.clone())
        }
    }

    
    /// 個数と値が組になった`Tuple`で解の値を返す
    pub fn value_tuple(&self) -> Vec<(usize, f64)> {
        match self {
            Root::Multiple(n, z) => vec![(*n, *z)],
            Root::Real(a1, a2) => {
                    let (a_l, a_s) =  Root::ord_z(a1, a2);
                    vec![a_l, a_s]
                },
            Root::Edge(n, z, r) => {
                    let mut vt = vec![(*n, *z)];
                    vt.append(&mut r.value_tuple());
                    vt
                },
        }
    }


    /// 解の分類を日本語のStringを出力
    pub fn kind_root(&self) -> String {
        match self {
            Root::Multiple(_,_) => "重根解".to_string(),
            Root::Real(a1, a2) => {
                let ((nl, _), (ns, _)) = Root::ord_z(a1, a2);
                format!("複号根解:{nl}/{ns}")
            },
            Root::Edge(e,_,r) => format!("エッジ解:K={}(擬{})", e, r.kind_root()),
        }
    }

    /// 解の値$ Z_i $の中に基準値以上の値を持つものが含まれるか確認
    ///
    /// 基準値よりも大きな値が含まれていれば`True`を返す．
    ///
    /// # 引数
    /// * `level` - 基準値．たとえばこれを値域の最大値$ (2 \alpha \delta + \beta) $で与えれば，解が有効な値であるか確認できる．
    pub fn greater(&self, level: &f64) -> bool {
        match self {
            Root::Multiple(_, z) => z > level,
            Root::Real((_, z1), (_, z2)) => (z1 > level) || (z2 > level),
            Root::Edge(_, z, r) => (z > level) || r.greater(level),
        }
    }

}

type Coefs = (f64, f64, f64);

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
    coef_eta: Coefs,
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

    // $ Z $の値を計算
    #[allow(dead_code)]
    fn calc_z(&self, z: &f64) -> f64 {
        let alpha = self.alpha();
        let beta = self.beta();
        
        if alpha == 0.0 && beta == 0.0 && z ==&0.0 {
            // θ=0の場合
            0.0
        } else {
            (z - beta) / ( 2.0 * alpha)
        }
    }

    // $ Z_i $ の値を[`Vec`]で返す
    #[allow(dead_code)]
    pub fn discriminant(&self) -> Vec<f64> {
        self.root
            .value()
            .par_iter()
            .map(|z| self.calc_z(z))
            .collect()
    }


    // $\eta(Z)$ の停留点の判定式に用いる値 $2 \alpha Z_i^{\ast} + \beta$の値を計算して，個数と値を組にした[`Tuple`]の[`Vec`]で返す．
    //
    // 解の種類によって個数が変わる可能性があるため注意
    #[allow(dead_code)]
    pub fn z_tuple(&self) -> Vec<(usize, f64)> {
        self.root
            .value_tuple()
            .par_iter()
            .map(|(n, z)| (*n, self.calc_z(z)) )
            .collect()
    }


    // $ 2 \alpha Z_i^{\ast} + \beta $の値を要素数3の[`Tuple`]で返す．
    //
    // 要素数3は解の値の種類が最大で3種類であるため．
    // 2種類以下で与えられる場合には，不足分は not a numberで補完する．
    #[allow(dead_code)]
    pub fn kind_roots(&self) -> (f64, f64, f64) {
        let roots = self.root.value_tuple();
        Self::extract_value_tuple(&roots)
    }


    // $ Z_i^{\ast} $の値を要素数3の[`Tuple`]で返す．
    //
    // 要素数3は解の値の種類が最大で3種類であるため．
    // 2種類以下で与えられる場合には，不足分は not a numberで補完する．
    #[allow(dead_code)]
    pub fn kind_z(&self) -> (f64, f64, f64) {
        let roots = self.z_tuple();
        Self::extract_value_tuple(&roots)
    }


    // 解の値を格納したタプルのベクトル`Vec<(usize,f64)>`から，値のみを抜き出して要素数3の[`Tuple`]で返す．
    //
    // 要素数3は解の値の種類が最大で3種類であるため．
    // 2種類以下で与えられる場合には，不足分は not a numberで補完する．
    fn extract_value_tuple(valt: &[(usize, f64)]) -> (f64, f64, f64) {
        let mut vals = valt.iter()
                           .map(|v| v.1)
                           .collect::<Vec<f64>>();
        let mut nans = vec![f64::NAN; 3 - vals.len()];
        vals.append(&mut nans);
        (vals[0], vals[1], vals[2])
    }
}



/// 有薗らによる多項式を分母に用いた確立不等式  
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
    /// 
    /// # 引数
    /// * `n` - 確率変数の個数
    /// * `sigma2` - 確率変数の分散の平均
    /// * `delta` - 確率変数の定義域の最大値
    /// 
    /// # 使用例
    /// ```
    /// # use hoeffding_ineq::polynominal2::PolynominalIneq;
    /// let polynominal = PolynominalIneq::new(4.0, 5.0).unwrap();
    /// ```
    pub fn new(n: usize, sigma2: f64, delta: f64) -> Result<Self, IneqError> {
        if sigma2 <= 0.0 {
            Err(IneqError{
                message: "sigma2 must be greater than zero.".to_string()
            })
        } else {
            Ok(PolynominalIneq{n, sigma2, delta})
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
    fn cal_coef_eta(&self, z: &f64) -> Coefs {
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


    // $ Z_i $から$ 2 \alpha Z_i + \beta $を計算
    fn root_val_z(z: &f64, coef_eta: &Coefs ) -> f64 {
        let (alpha, beta, _) = coef_eta;
        2.0 * alpha * z + beta
    }


    // 重根解の計算
    // 
    // # 引数
    // * `n` - 確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn root_multiple(n: &usize, z: &f64, coef_eta: &Coefs) -> Root {
        let val = Self::root_val_z(z, coef_eta);
        Root::Multiple(*n, val)
    }


    // 複号根解の計算
    //
    // # 引数
    // * `n` - 確率変数の個数
    // * `m` - 解のうち大きな値である解$ Z_{(2)} $の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    //
    // # 返り値について
    // 複号根解が存在しない場合[`Option::None`]を返す．
    // 複号根解が存在する場合は[`Option::Some`]を返す．
    // ただし，解$ Z $の値が適切な値域に入っているかは確認していないため注意．
    fn set_root_real(n: &usize, m: &usize, z: &f64, coef_eta: &Coefs) -> Option<Root> {
        let n_f = *n as f64;
        let m_f = *m as f64;
        let n_m_f = (*n - *m) as f64;
        let (alpha, beta, gamma) = coef_eta;
        let sum_val = Self::root_val_z(z, coef_eta) * n_f;
        let discri = sum_val.powi(2) - (4.0 * n_m_f * m_f * (4.0 * alpha * gamma - (beta.powi(2)) ) );

        if discri.is_sign_negative() {
            None
        } else {
            let z1 = (sum_val - discri.sqrt()) / (2.0 * n_m_f);
            let z2 = (sum_val + discri.sqrt()) / (2.0 * m_f);
            Some( Root::Real((n-m, z1), (*m, z2)) )
        }
    }


    // 複号根解の計算
    // 
    // # 引数
    // * `n` - 確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    //
    // # 返り値について
    // 実数解を与えるもののみを抽出．
    // ただし，解$ Z_i $の値が適切な値域に入っているかは確認していないため注意．
    fn root_real(n: &usize, z: &f64, coef_eta: &Coefs ) -> Vec<Root> {
        (1..*n).collect::<Vec<usize>>()
               .par_iter()
               .filter_map(|m| Self::set_root_real(n, m, z, coef_eta) )
               .collect()
    }


    // 複号根解の中で$ \prod \eta(Z_i) $の最小値を与える解を抽出
    //
    // 解が存在しない場合があり，その場合は[`Option::None`]を返す．
    //
    // # 引数
    // * `n` - 確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    //
    // # 注意
    // 解$ Z_i $の値が適切な値域に入っているかは確認していないため注意．
    fn root_real_min(n: &usize, z: &f64, coef_eta: &Coefs) -> Option<Root> {
        let roots = Self::root_real(n, z, coef_eta);
        if roots.len() == 0 {
            None
        } else {
            Some( Self::argmin_root(coef_eta, &roots) )
        }
    }


    // 複号根解について，全ての解の値域が適切かを確認し，値域が適切なものを抽出する．
    //
    // # 引数
    // * `n` - 確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `delta` - $ Z_i $の値域の上限
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    //
    // # 返り値
    // `(out_of_range, roots)`  
    // * `out of range` - 全ての解が値域の上限$ \delta $を満たしているか．満たしていない解が存在する場合は`false`．
    // * `roots` - 条件を満たした解
    fn root_real_in_range(&self, n: &usize, z: &f64, coef_eta: &Coefs ) -> Vec<Root> {
        let roots = Self::root_real(n, z, coef_eta);

        let level = self.root_delta(coef_eta);

        if roots.par_iter().any(|r| r.greater(&level)) {
            let roots_in_range = roots.into_par_iter()
                                      .filter_map(|r| 
                                            if r.greater(&level) { 
                                                None
                                            } else {
                                                Some(r)
                                            }
                                     ).collect::<Vec<Root>>();
            
            roots_in_range
        } else {
            roots
        }        
    }


    // 解の値が確率変数の上限値を満たす複号根解の中で$ \prod \eta(Z_i) $の最小値を与える解を抽出
    //
    // 解が存在しない場合があり，その場合は[`Option::None`]を返す．
    //
    // # 引数
    // * `n` - 確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn root_real_min_in_range(&self, n: &usize, z: &f64, coef_eta: &Coefs) -> Option<Root> {
        let roots = self.root_real_in_range(n, z, coef_eta);
        if roots.len() == 0{
            // 値域の条件を満たす解なし
            None
        } else {
            Some( Self::argmin_root(coef_eta, &roots) )
        }
    }


    // 解$ 2 \alpha Z_i + \beta $の上限値$ 2 \alpha \delta + \beta $
    //
    // # 引数
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn root_delta(&self, coef_eta: &Coefs ) -> f64 {
        Self::root_val_z(&self.delta, coef_eta)
    }



    // $ K $個の変数を$ \delta $で固定した場合のエッジ解の計算
    //
    // # 引数
    // * `k` - 値を$ Z_i = \delta $で固定する確率変数の個数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn set_root_edge(&self, k: &usize, z: &f64, coef_eta: &Coefs) -> Option<Root> {
        let n_f = self.n as f64;
        let level = self.root_delta(coef_eta);

        let new_n = self.n - k;
        let new_z = (n_f * z - (*k as f64) * self.delta) / (new_n as f64);

        let min_statio = self.calc_stationary_min_in_range(&new_n, &new_z, coef_eta);

        // 解に$ Z_i = \delta $を含むか確認
        if min_statio.value().iter().any(|r| *r == level) {
            None
        } else {
            Some( Root::Edge(*k, level, Box::new(min_statio)) )
        }
    }


    // エッジ解の計算
    //
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn root_edge(&self, z: &f64, coef_eta: &Coefs) -> Root {
        let (alpha, beta, _) = coef_eta;
        let n_f = self.n as f64;
        let k_max = ( (n_f * ( 2.0 * alpha * (*z) + beta)) / (2.0 * alpha * self.delta + beta) ).floor() as usize;

        let roots_edge = (1..=k_max).into_iter()
                                    .filter_map(|k| self.set_root_edge(&k, z, coef_eta))
                                    .collect::<Vec<Root>>();
        Self::argmin_root(coef_eta, &roots_edge)
    }


    // 評価値$ \prob \eta(Z_i) $を最小にする停留点（重根解 or エッジ解）を返す
    // 
    // 内部計算用の関数であり，引数に$ \eta(Z_i) $の係数$ \alpha, \beta, \gamma $を取る
    // この関数はエッジ解以外の解で最小の値を返す．
    // 解の範囲のチェック無し
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn calc_stationary_min(n: &usize, z: &f64, coef_eta: &Coefs) -> Root {
        let mult = Self::root_multiple(n, z, coef_eta);
        let real_opt = Self::root_real_min(n, z, coef_eta);
        if let Some(real) = real_opt {
            Self::argmin_root(coef_eta, &[mult, real])
        } else {
            // 複号根解が存在しない場合
            mult
        }
    }


    // 評価値$ \prob \eta(Z_i) $を最小にする$ \delta $以下の停留点（重根解 or エッジ解）を返す
    // 
    // 内部計算用の関数であり，引数に$ \eta(Z_i) $の係数$ \alpha, \beta, \gamma $を取る
    // この関数はエッジ解以外の解で最小の値を返す．
    // 解の範囲のチェックあり
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    fn calc_stationary_min_in_range(&self, n: &usize, z: &f64, coef_eta: &Coefs) -> Root {
        let mult = Self::root_multiple(n, z, coef_eta);
        let real_opt = self.root_real_min_in_range(n, z, coef_eta);
        if let Some(real) = real_opt {
            Self::argmin_root(coef_eta, &[mult, real])
        } else {
            // 複号根解が存在しない場合
            mult
        }
    }


    // 評価値$ \prob \eta(Z_i) $を最小とする解を探索
    //
    // 内部計算用の関数であり，引数に$ \eta(Z_i) $の係数$ \alpha, \beta, \gamma $を取る
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    #[allow(dead_code)]
    fn calc_optimal_root(&self, z: &f64, coef_eta: &Coefs) -> Root {
        if *z == self.delta {
            return Root::Multiple(self.n, self.root_delta(coef_eta));
        }

        let min_root = Self::calc_stationary_min(&self.n, z, coef_eta);
        if min_root.greater(&self.root_delta(coef_eta)) {
            // エッジ解が必要な場合
            let min_inrange = self.calc_stationary_min_in_range(&self.n, z, coef_eta);
            let min_edge = self.root_edge(z, coef_eta);
            Self::argmin_root(coef_eta, &[min_inrange, min_edge])
        } else {
            min_root
        }
    }


    /// 与えられた変数 $z$ に対して，確立不等式の分母となる関数 $\eta(Z)$ の停留点を与える解
    /// 
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn optimal_root(&self, z: &f64) -> Root {
        let coef_eta = self.cal_coef_eta(z);
        self.calc_optimal_root(z, &coef_eta)
    }


    // $ e^{Z_i \theta / n) $を上から押さえつける関数$ \eta(Z_i) $の計算
    //
    // # 引数
    // * `root_val` - $ 2 \alpha Z_i + \beta $の値 
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    // * `root` - $ prod ta(Z_i) $の最小値を与える解の候補
    fn calc_eta(root_val: &f64, coef_eta: &Coefs) -> f64 {
        let (alpha, beta, gamma) = coef_eta;
        (   root_val.powi(2)
            +
            4.0 * alpha * gamma - beta.powi(2)
        ) / (4.0 * alpha)
    }


    // 確立不等式の分母となる，停留点での $\eta(Z_i)$ の値の積
    // 内部計算用のため，すべての変数を引数に取ります．
    // 利用する際は適宜ラッパー関数を利用してください．
    // 
    // # 引数
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    // * `root` - $ prod ta(Z_i) $の最小値を与える解の候補
    fn calc_pi_eta(coef_eta: &Coefs, root: &Root) -> f64 {
        match root {
            Root::Multiple(n, rv) => {
                Self::calc_eta(rv, coef_eta).powf(*n as f64)
            },
            Root::Real((n1, rv1), (n2, rv2)) => {
                Self::calc_eta(rv1, coef_eta).powf(*n1 as f64)
                *
                Self::calc_eta(rv2, coef_eta).powf(*n2 as f64)
            },
            Root::Edge(n, rv, root_rest) => {
                Self::calc_pi_eta(coef_eta, root_rest)
                *
                Self::calc_eta(rv, coef_eta).powf(*n as f64)
            }
        }
    }


    // $ \prod \eta(Z_i) $の最小値を与える[`Root`]を選んで返す
    //
    // # 注意
    // 引数に与えるRootの[`Vec`]が空の場合はpanicになります．
    //
    // # 引数
    // * `coef_eta` - 確率不等式の分母となる関数 $\eta(Z)$ の係数
    //      0. `alpha` - $\eta(Z)$ の2次の項の係数
    //      1. `beta` - $\eta(Z)$ の1次の項の係数
    //      2. `gamma` - $\eta(Z)$ の定数項
    // * `roots` - $ prod ta(Z_i) $の最小値を与える解の候補
    fn argmin_root(coef_eta: &Coefs, roots: &[Root]) -> Root {
        let vals_prod_eta = roots.par_iter()
                                 .map(|r| Self::calc_pi_eta(coef_eta, r))
                                 .collect::<Vec<f64>>();
        let zip_r_v = std::iter::zip(roots, vals_prod_eta);
        let min_r_v = zip_r_v.reduce(|(r_acc, v_acc), (r_i, v_i)|
                                if v_acc.lt(&v_i) {
                                    (r_acc, v_acc)
                                } else {
                                    (r_i, v_i)
                                }
                            ).unwrap();
        min_r_v.0.clone()
    }


    // 上側確率の上界を計算
    // 内部計算用のため，すべての変数を引数に取ります．
    // 利用する際は適宜ラッパー関数を利用してください．
    // 
    // # 引数
    // * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値 
    fn calculate(&self, z: &f64) -> (f64, Root) {
        if z == &0.0 {
            // z = 0 ならθ=0であり分子は1.0，さらにα=0,β=0,γ=1よりη(Z)=1が成り立つ
            // よって必ず1.0となる
            return (1.0, Root::Multiple(self.n, 0.0))
        }

        let coef_eta = self.cal_coef_eta(z);
        let root = self.calc_optimal_root(z, &coef_eta);
        let pi_eta = Self::calc_pi_eta(&coef_eta, &root);
        let theta = self.theta(z);
        let d2ps2 = self.delta.powi(2) + self.sigma2;

        let prob = ((self.delta.powi(2) / d2ps2 * ((- theta / (self.n as f64)* self.sigma2 / self.delta).exp())
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
        if *z < self.min_z() {
            Err(IneqError {
                message: format!("`z` must be greater {}", self.min_z())
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
        let pi_eta = if z == &0.0 {
            // z=0ならθ=0,α=0,β=0,γ=1よりη=1
            // よってΠη(z)=1
            1.0
        } else {
            Self::calc_pi_eta(&coef_eta, &root)
        };
        let theta = self.theta(z);
        Ok(VerbosePolynominalIneq{z: *z, pr, root, coef_eta, pi_eta, theta})
    }


    // 引数の有効な範囲内における上側確率の上界を，詳細情報とともにまとめて計算  
    fn overview_with_verbose(&self) -> Result<Vec<VerbosePolynominalIneq>, IneqError> {
        let vec_z = self.variable_for_overview();
        vec_z.par_iter()
             .map(|z| self.pr_with_verbose(z))
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
            sheet_writer.append_row(xlsx::row!["z", "提案手法","2αZ1+β", "2αZ2+β", "2αZ3+β", "Z1", "Z2", "Z3", "根の種類"])?;
            for ver in overview.iter() {
                let (z1, z2, z3) = ver.kind_z();
                let (root1, root2, root3) = ver.kind_roots();
                // NANの削除
                let (z2_cell, root2_cell) = if f64::is_nan(root2) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z2.to_cell_value(), root2.to_cell_value())
                        };
                let (z3_cell, root3_cell) = if f64::is_nan(root3) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z3.to_cell_value(), root3.to_cell_value())
                        };
                sheet_writer.append_row(xlsx::row![
                        ver.z.to_cell_value(), 
                        ver.pr.to_cell_value(), 
                        root1.to_cell_value(),
                        root2_cell,
                        root3_cell,
                        z1.to_cell_value(),
                        z2_cell,
                        z3_cell,
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
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Z1", "Z2", "Z3", "2αZ1+β", "2αZ2+β", "2αZ3+β", "種類", "theta", "alpha", "beta", "gamma", "Pi eta(Z_i)","Hoeffdingの分母","相対改善比率"])?;
            for ver in overview.iter() {
                let (z1, z2, z3) = ver.kind_z();
                let (root1, root2, root3) = ver.kind_roots();
                // NANの削除
                let (z2_cell, root2_cell) = if f64::is_nan(root2) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z2.to_cell_value(), root2.to_cell_value())
                        };
                let (z3_cell, root3_cell) = if f64::is_nan(root3) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z3.to_cell_value(), root3.to_cell_value())
                        };
                
                let hoef_ezt = (ver.theta * ver.z).exp();
                let improvement = (ver.pi_eta - hoef_ezt) / ver.pi_eta;
                sheet_writer.append_row(xlsx::row![
                    ver.z.to_cell_value(), 
                    ver.pr.to_cell_value(),
                    z1.to_cell_value(),
                    z2_cell,
                    z3_cell,
                    root1.to_cell_value(),
                    root2_cell,
                    root3_cell,
                    ver.root.kind_root(),
                    ver.theta.to_cell_value(),
                    ver.alpha().to_cell_value(), 
                    ver.beta().to_cell_value(), 
                    ver.gamma().to_cell_value(),
                    ver.pi_eta.to_cell_value(),
                    hoef_ezt.to_cell_value(),
                    improvement.to_cell_value()
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
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Hoeffding", "相対改善比率","2αZ1+β", "2αZ2+β", "2αZ3+β", "根の種類"])?;
            for ver in overview.iter() {
                let (root1, root2, root3) = ver.kind_roots();
                let pr_hoeff = hoeff.pr(&ver.z);
                // NANの削除
                let root2_cell = if f64::is_nan(root2) {
                            ().to_cell_value()
                        } else {
                            root2.to_cell_value()
                        };
                let root3_cell = if f64::is_nan(root3) {
                            ().to_cell_value()
                        } else {
                            root3.to_cell_value()
                        };
                let ratio = (pr_hoeff - ver.pr) / pr_hoeff;

                sheet_writer.append_row(xlsx::row![
                        ver.z.to_cell_value(), 
                        ver.pr.to_cell_value(), 
                        pr_hoeff.to_cell_value(),
                        ratio.to_cell_value(),
                        root1.to_cell_value(),
                        root2_cell,
                        root3_cell,
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
            sheet_writer.append_row(xlsx::row!["z", "提案手法", "Z1", "Z2", "Z3", "2αZ1+β", "2αZ2+β", "2αZ3+β", "種類", "theta", "alpha", "beta", "gamma", "Pi eta(Z_i)"])?;
            for ver in overview.iter() {
                let (z1, z2, z3) = ver.kind_z();
                let (root1, root2, root3) = ver.kind_roots();
                // NANの削除
                let (z2_cell, root2_cell) = if f64::is_nan(root2) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z2.to_cell_value(), root2.to_cell_value())
                        };
                let (z3_cell, root3_cell) = if f64::is_nan(root3) {
                            (().to_cell_value(), ().to_cell_value())
                        } else {
                            (z3.to_cell_value(), root3.to_cell_value())
                        };
                sheet_writer.append_row(xlsx::row![
                    ver.z.to_cell_value(), 
                    ver.pr.to_cell_value(),
                    z1.to_cell_value(),
                    z2_cell,
                    z3_cell,
                    root1.to_cell_value(),
                    root2_cell,
                    root3_cell,
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


    // `Root`とその評価値の`Root`に依存する部分($\prod_{i=1}^{n} \eta(Z_i^{\ast})$)の組を作成
    //
    // # 返り値
    // * `kind, vals, eval`
    //  * `kind` - 解の種類
    //  * `vals` - 解の個数と値の組．要素数は3で，足りない分は`(0, f64::NAN)`で補完
    //  * `eval` - 評価値の`Root`に依存する部分($\prod_{i=1}^{n} \eta(Z_i^{\ast})$) 
    fn root_and_pi_eta(coef_eta: &Coefs, root: &Root) -> (Root, Vec<(usize, f64)>, f64) {
        let mut vals = root.value_tuple();
        let mut nans = vec![(0, f64::NAN); 3-vals.len()];
        vals.append(&mut nans);
        let eval = Self::calc_pi_eta(coef_eta, root);
        (root.clone(), vals, eval)
    }


    /// $z$の値に対して，候補となる各根での値を取得
    ///
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn info_each_root(&self, z: f64) -> Result<Vec<(Root, Vec<(usize, f64)>, f64)>, IneqError> {
        let coef_eta = self.cal_coef_eta(&z);
        // 重根解
        let mult = Self::root_multiple(&self.n, &z, &coef_eta);
        let mut each_root = vec![Self::root_and_pi_eta(&coef_eta, &mult)];
        // 複号根解
        let mut root_real = Self::root_real(&self.n, &z, &coef_eta)
                                 .par_iter()
                                 .rev()
                                 .map(|r| Self::root_and_pi_eta(&coef_eta, r))
                                 .collect();
        each_root.append(&mut root_real);
        // エッジ解
        let min_root = Self::calc_stationary_min(&self.n, &z, &coef_eta);
        if min_root.greater(&self.root_delta(&coef_eta)) {
            // エッジ解が必要な場合
            let min_edge = self.root_edge(&z, &coef_eta);
            each_root.push(Self::root_and_pi_eta(&coef_eta, &min_edge));
        };
        Ok(each_root)
    }


    /// $z$の値に対して，候補となる各根での値を取得し，Excelファイルに保存
    ///
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn info_each_root_to_excel(&self, z: f64, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {

        // Excelへの書き込み

        // 根の情報
        let mut wb = xlsx::Workbook::create(xlsx_path);
        let mut sheet_1 = wb.create_sheet("Roots");
        
        sheet_1.add_column(xlsx::Column { width: 28.0 });
        for _i in 1..=9 {
            sheet_1.add_column(xlsx::Column { width: 11.0 });
        }
        for _i in 11..=16 {
            sheet_1.add_column(xlsx::Column { width: 15.5 });
        }

        let info_roots = self.info_each_root(z)?;
        let verbose = self.pr_with_verbose(&z)?;
        wb.write_sheet(&mut sheet_1, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["根の種類", "Z1の個数", "2αZ1+β", "Z1", "Z2の個数", "2αZ2+β", "Z2", "Z=δの個数", "2αδ+β", "Π η(Z*)", "エッジ解:Z1の個数", "エッジ解:2αZ1+β", "エッジ解:Z1", "エッジ解:Z2の個数", "エッジ解:2αZ2+β", "エッジ解:Z2"])?;
            for ier in info_roots {
                let (root, vals, eval) = ier;
                match root {
                    Root::Edge(_, _, _) => {
                        let z2_n_cell;
                        let z2_val_cell;
                        let z2_cell;
                        if f64::is_nan(vals[2].1) {
                            z2_n_cell = ().to_cell_value();
                            z2_val_cell = ().to_cell_value();
                            z2_cell = ().to_cell_value();
                        } else {
                            z2_n_cell = (vals[2].0 as f64).to_cell_value();
                            z2_val_cell = vals[2].1.to_cell_value();
                            z2_cell = verbose.calc_z(&vals[2].1).to_cell_value();
                        }
                        // エッジ解の場合には順序を入れ替える
                        sheet_writer.append_row(xlsx::row![root.kind_root().to_cell_value(),
                                                           (),
                                                           (),
                                                           (),
                                                           (),
                                                           (),
                                                           (),
                                                           (vals[0].0 as f64).to_cell_value(),
                                                           vals[0].1.to_cell_value(),
                                                           eval.to_cell_value(),
                                                           (vals[1].0 as f64).to_cell_value(),
                                                           vals[1].1.to_cell_value(),
                                                           verbose.calc_z(&vals[1].1).to_cell_value(),
                                                           z2_n_cell,
                                                           z2_val_cell,
                                                           z2_cell])?
                        },
                    Root::Multiple(_, _) => {
                        sheet_writer.append_row(xlsx::row![root.kind_root().to_cell_value(),
                                                           (vals[0].0 as f64).to_cell_value(),
                                                           vals[0].1.to_cell_value(),
                                                           verbose.calc_z(&vals[0].1).to_cell_value(),
                                                           (),
                                                           (),
                                                           (),
                                                           (),
                                                           (),
                                                           eval.to_cell_value()])?
                        },
                    Root::Real(_, _) => {
                        sheet_writer.append_row(xlsx::row![root.kind_root().to_cell_value(),
                                                           (vals[0].0 as f64).to_cell_value(),
                                                           vals[0].1.to_cell_value(),
                                                           verbose.calc_z(&vals[0].1).to_cell_value(),
                                                           (vals[1].0 as f64).to_cell_value(),
                                                           vals[1].1.to_cell_value(),
                                                           verbose.calc_z(&vals[1].1).to_cell_value(),
                                                           (),
                                                           (),
                                                           eval.to_cell_value()])?
                        },
                };
            }
            // Hoeffdingの確率不等式による分母を記載
            let hoef_ezt = (z * verbose.theta).exp();
            sheet_writer.append_row(xlsx::row!["Hoeffdingの確率不等式", (), (), (), (), (), (), (), (), hoef_ezt.to_cell_value()])?;
            Ok(())
        })?;

        // シミュレーション情報
        let mut sheet_2 = wb.create_sheet("Parameters");
        
        let params = self.param_to_tuple();
        wb.write_sheet(&mut sheet_2, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["Parameter", "Value"])?;
            for (p, v) in params.iter() {
                sheet_writer.append_row(xlsx::row![p.to_cell_value(), v.to_cell_value()])?;
            }
            sheet_writer.append_row(xlsx::row!["z".to_cell_value(), z.to_cell_value()])?;
            sheet_writer.append_row(xlsx::row!["上界".to_cell_value(), verbose.pr.to_cell_value()])?;
            sheet_writer.append_row(xlsx::row!["theta".to_cell_value(), verbose.theta.to_cell_value()])?;

            Ok(())
        })?;
    
        Ok(())
    }
}

impl UpperBound for PolynominalIneq {
    fn pr(&self, z: &f64) -> f64 {
        if *z < self.min_z() {
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
