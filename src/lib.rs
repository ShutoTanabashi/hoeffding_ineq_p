//! 上側確率の上界を与えるHoeffdingの確率不等式に関する計算

pub mod hoeffding;
pub mod polynominal;
pub mod polynominal2;
pub mod compare_theoretical;

use std::{self, fmt};
extern crate rayon;
use rayon::prelude::*;
extern crate simple_excel_writer;
use simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};


/// 確率不等式に関するエラー
#[derive(Debug, Clone)]
pub struct IneqError {
    pub message: String,
}

impl fmt::Display for IneqError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for IneqError {
    fn description(&self) -> &str {
        &self.message
    }
}


/// 上側確率の上界を計算するトレイト
pub trait UpperBound: Sync {
    /// 上側確率の上界を計算  
    /// 
    /// # 注意
    /// 本来は計算範囲外である $z < 0$ や $z > \delta$ の場合も適宜値を返す．
    /// 
    /// # 引数  
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    fn pr(&self, z: &f64) -> f64;

    /// 引数zの値域の上限
    fn max_z(&self) -> f64;

    /// 引数zの値域の下限
    fn min_z(&self) -> f64;

    /// 各種パラメータをタプルで出力する
    /// 全てのパラメータが浮動小数点数で出力される点に注意．
    fn param_to_tuple(&self) -> Vec<(String, f64)>;

    /// 上側確率の上界を計算し，(z, Pr(z))のタプルを返す
    /// 
    /// # 引数  
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    fn pr_tuple(&self, z: &f64) -> (f64, f64) {
        let pr = self.pr(z);
        (z.clone(), pr)
    }


    /// overview や variable_for_overview の返り値の要素数  
    /// すなわち，overviewにてどれほど細かく計算するかを示す．
    const NUM_POINT_OVERVIEW: u64 = 500;

    /// 複数個の変数に対して上側確率の上界を計算
    /// 
    /// # 引数
    /// * `zs` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値の組
    fn map_pr(&self, zs: &[f64]) -> Vec<f64> {
        zs.par_iter()
          .map(|z| self.pr(z))
          .collect::<Vec<f64>>()
    }

    
    /// 引数の有効な範囲について，一定間隔で確率変数を生成する
    fn variable_for_overview(&self) -> Vec<f64> {
        let gap_point = (self.max_z() - self.min_z()) / (Self::NUM_POINT_OVERVIEW as f64);
        (0..Self::NUM_POINT_OVERVIEW).map(|i| (i as f64) * gap_point)
                               .collect::<Vec<f64>>()
    }


    /// 引数の有効な範囲内における上側確率の上界をまとめて計算
    /// 返り値のタプルは(z, Pr(z))を意味する
    fn overview(&self) -> Vec<(f64, f64)> {
        let vec_z = self.variable_for_overview();
        vec_z.par_iter()
             .map(|z| self.pr_tuple(z))
             .collect::<Vec<(f64, f64)>>()
    }

    /// 引数の有効な範囲内における上側確率の上界をExcelファイルで保存
    /// 
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    fn overview_to_excel(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let overview = self.overview();
        let params = self.param_to_tuple();

        let mut wb = xlsx::Workbook::create(xlsx_path);
        let mut sheet_1 = wb.create_sheet("Calculated Value");

        wb.write_sheet(&mut sheet_1, |sheet_writer| {
            sheet_writer.append_row(xlsx::row!["z", "Pr(Z>=z)"])?;
            for (z, prz) in overview.iter() {
               sheet_writer.append_row(xlsx::row![z.to_cell_value(), prz.to_cell_value()])?;
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

/// 上側確率の上界をHoeffdingの確率不等式と比較可能
pub trait CompareHoeffding: UpperBound {
    /// 同条件のHoeffdingの確率不等式を生成する
    fn same_condition_hoeffding(&self) -> Result<hoeffding::HoeffdingIneq, IneqError>;

    /// 変数 z に対して，SelfとHoeffdingの確率不等式の両手法での上側確率の上界を比較する
    ///
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    ///
    /// # 返り値
    /// * `(z, self_pr, hoef_pr)` - 上側確率の上界
    ///     * `self_pr` - Selfで計算された上側確率の上界
    ///     * `hoef_pr` - Hoeffdingの確率不等式で計算された上側確率の上界
    fn compare_hoeffding(&self, z: &f64) -> Result<(f64, f64, f64), IneqError> {
        let hoef = self.same_condition_hoeffding()?;
        let hoef_pr = hoef.pr(z);
        let self_pr = self.pr(z);
        Ok((z.clone(), self_pr, hoef_pr))
    }
}

