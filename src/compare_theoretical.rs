//! 理論値との比較
extern crate statrs;
use statrs::distribution::ChiSquared;
use statrs::distribution::ContinuousCDF;
use statrs::StatsError;
extern crate rayon;
use rayon::prelude::*;
extern crate simple_excel_writer;
use simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};

use super::{UpperBound, CompareHoeffding, IneqError};
use super::polynominal::PolynominalIneq;
use super::hoeffding::HoeffdingIneq;


/// 有薗ら提案の多項式に基づくHoeffding型確率不等式およびHoeffdingの確率不等式に対して，理論値を用いた比較が可能
pub trait ComparePolynominal: UpperBound {
    /// 乱数の個数を取得
    fn get_n(&self) -> usize;
    /// 乱数の分散の代表値を取得(`E[σ_i]`)
    fn get_sigma2(&self) -> f64;
    /// 乱数の最大値の代表値を取得(`max δ_i`)
    fn get_delta(&self) -> f64;


    /// 確率変数の条件に合致する2個の確率不等式を取得
    ///
    /// # 返り値
    /// * `(poly_ineq, hoef_ineq)`
    ///     * `poly_ineq` - 有薗ら提案の多項式に基づくHoeffding型確率不等式
    ///     * `hoef_ineq` - Hoeffdingの確率不等式
    fn same_condition_ineqs(&self) -> Result<(PolynominalIneq, HoeffdingIneq), IneqError> {
        let poly_ineq = self.same_condition_polynominal()?;
        let hoef_ineq = poly_ineq.same_condition_hoeffding()?;
        Ok((poly_ineq, hoef_ineq))
    }


    /// 確率変数の条件に合致する，有薗ら提案の多項式に基づくHoeffding型確率不等式を取得
    fn same_condition_polynominal(&self) -> Result<PolynominalIneq, IneqError> {
        PolynominalIneq::new(self.get_n(), self.get_sigma2(), self.get_delta())
    }


    /// 変数 z に対して，理論値と有薗らによる多項式に基づく確率不等式，Hoeffdingの確率不等式の3種類の手法での上側確率の上界を比較する
    ///
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    ///
    /// # 返り値
    /// * `(z, poly_pr, hoef_pr, theo_pr)` - 上側確率の上界
    ///     * `theo_pr` - 上側確率の理論値
    ///     * `poly_pr` - 有薗らによる多項式に基づくHoeffding型確率不等式で計算された上側確率の上界
    ///     * `hoef_pr` - Hoeffdingの確率不等式で計算された上側確率の上界
    fn compare(&self, z: &f64) -> Result<(f64, f64, f64, f64), IneqError> {
        let (poly_ineq, hoef_ineq) = self.same_condition_ineqs()?;
        let theo_pr = self.pr(z);
        let poly_pr = poly_ineq.pr(z);
        let hoef_pr = hoef_ineq.pr(z);
        Ok((z.clone(), theo_pr, poly_pr, hoef_pr))
    }


    /// 引数の有効な範囲内における上側確率の上界をまとめて計算
    ///
    /// 範囲は有薗らによるhoeffding型確率不等式(`polynominal`)に合わせる
    ///
    /// # 返り値
    /// * `Vec<(z, poly_pr, hoef_pr, theo_pr)>` - 上側確率の上界
    ///     * `theo_pr` - 上側確率の理論値
    ///     * `poly_pr` - 有薗らによる多項式に基づくHoeffding型確率不等式で計算された上側確率の上界
    ///     * `hoef_pr` - Hoeffdingの確率不等式で計算された上側確率の上界
    fn compare_overview(&self) -> Result<Vec<(f64, f64, f64, f64)>, IneqError> {
        let (poly_ineq, hoef_ineq) = self.same_condition_ineqs()?;
        let vec_z = poly_ineq.variable_for_overview();
        let prs = vec_z.par_iter()
                       .map(|z| (z.clone(), self.pr(z), poly_ineq.pr(z), hoef_ineq.pr(z)) )
                       .collect::<Vec<(f64, f64, f64, f64)>>();
        Ok(prs)
    }



    /// 引数の有効な範囲内における上側確率の上界をExcelファイルで保存
    ///
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    fn compare_overview_to_excel(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let overview = self.compare_overview()?;
        let params = self.param_to_tuple();

        let mut wb = xlsx::Workbook::create(xlsx_path);
        let mut sheet_1 = wb.create_sheet("Calculated Value");

        wb.write_sheet(&mut sheet_1, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "理論値", "提案手法", "Hoeffding", "相対改善比率"])?;

            for (z, theo_pr, poly_pr, hoef_pr) in overview.iter() {
                let improvement = (hoef_pr - poly_pr) / hoef_pr;
                sheet_writer.append_row(xlsx::row![z.to_cell_value(), theo_pr.to_cell_value(), poly_pr.to_cell_value(), hoef_pr.to_cell_value(), improvement.to_cell_value()])?;
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

        Ok(())
    }
}



/// カイ二乗分布の理論値との比較を実行
pub fn test_with_chi2(dfs: &[f64]) {
    println!("start compare with chi_2 for {:?}", dfs);
    let chi2 = TheoreticalChiSquared::new(dfs).unwrap();
    // let vals = [0.0 , 0.5, 1.0, 1.5, 2.0];
    // let sf = vals.into_iter()
    //              .map(|v| chi2.compare(&v))
    //              .collect::<Result<Vec<(f64,f64,f64,f64)>,IneqError>>()
    //              .unwrap();
    let sf = chi2.compare_overview().unwrap();
    println!("{:?}", sf);

    println!("{:?}", chi2.param_to_tuple());
}


/// カイ二乗分布に従う変数の和の理論値
///
/// 自由度`d_i`のカイ二乗分布に従う変数`x_i`に対して，`z_i = -x_i + d_i`で変数変換を行い，`z_i`の平均`Z`が基準値`z`を上回る確率を求める．
/// `Pr{Z>=z} = `Pr{Σz_i < Σ(d_i - z)}`
pub struct TheoreticalChiSquared {
    pub dfs: Vec<f64>,
    chi2: ChiSquared,
}

impl TheoreticalChiSquared {
    /// 各変数の自由度から理論値計算用のインスタンスを生成
    ///
    /// # 引数
    /// * `dfs` - 確率変数が従うそれぞれのカイ二乗分布の自由度
    pub fn new(dfs: &[f64]) -> Result<Self,StatsError> {
        let sum = dfs.iter()
                     .sum();
        let chi2 = ChiSquared::new(sum)?;
        let dfs_vec = dfs.to_vec();
        Ok(TheoreticalChiSquared { dfs: dfs_vec, chi2 })
    }

    /// 任意のzに対して上側確立を計算
    ///
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    pub fn sf(&self, z: f64) -> f64 {
        let z_t = self.dfs
                      .iter()
                      .map(|df| df - z)
                      .sum();
        self.chi2
            .cdf(z_t)
    }
}


impl UpperBound for TheoreticalChiSquared {
    fn pr(&self, z: &f64) -> f64 {
        self.sf(*z)
    }

    fn max_z(&self) -> f64 {
        self.dfs
            .iter()
            .copied()
            .fold(f64::NAN, f64::max)
    }

    fn min_z(&self) -> f64 {
        f64::NEG_INFINITY
    }

    fn param_to_tuple(&self) -> Vec<(String, f64)> {
        let mut info_chi2 = vec![("df_theoretical".to_owned(), self.chi2.freedom()), 
                                 ("n".to_owned(), self.dfs.len() as f64)];
        let mut vec_df = self.dfs
                             .iter()
                             .enumerate()
                             .map(|(i, df)| (format!("df_{}", i), *df))
                             .collect::<Vec<(String, f64)>>();
        info_chi2.append(&mut vec_df);
        info_chi2
        
    }
}


impl ComparePolynominal for TheoreticalChiSquared {
    fn get_n(&self) -> usize {
        self.dfs.len()
    }

    fn get_sigma2(&self) -> f64 {
        let n = self.get_n() as f64;
        let cum = self.dfs
                      .iter()
                      .fold(0.0, |sum, x| sum + (2.0 * x) );
        cum / n
    }

    fn get_delta(&self) -> f64 {
        self.dfs
            .iter()
            .copied()
            .fold(0.0, f64::max)
    }
}

