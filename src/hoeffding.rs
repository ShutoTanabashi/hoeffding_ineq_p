//! Hoeffdingの確率不等式の計算

use super::{IneqError, UpperBound};


/// Hoeffdingの確率不等式
#[derive(Debug, Clone)]
pub struct HoeffdingIneq {
    // 確率変数の個数
    n: usize,
    // 確率変数の分散
    sigma2: f64,
    // 確率変数の最大値
    delta: f64,
}


impl HoeffdingIneq {
    /// 各種引数からHoeffdingの確率不等式を初期化
    ///
    /// # 引数
    /// * `n` - 確率変数の個数
    /// * `sigma2` - 確率変数の分散の平均
    /// * `delta` - 確率変数の定義域の最大値
    /// 
    /// # 使用例
    /// ```
    /// # use hoeffding_ineq::hoeffding::HoeffdingIneq;
    /// let hoeffding = HoeffdingIneq::new(5, 4.0, 5.0).unwrap();
    /// ```
    pub fn new(n: usize, sigma2: f64, delta: f64) -> Result<Self, IneqError> {
        if sigma2 <= 0.0 {
            Err(IneqError{
                message: "sigma2 must be greater than zero.".to_string()
            })
        } else {
            Ok(HoeffdingIneq{n, sigma2, delta})
        }
    }

    /// Hoeffdingの確立不等式に基づいて上側確率の上界を計算する
    /// ただし，引数の定義域の確認はなし．
    /// 定義域に応じて計算処理を変えるラッパー関数を適宜利用すること．
    /// 
    /// # 引数
    /// * `z` - $Pr\{\bar{Z} \leq z\}$ すなわち計算したい確率変数の値  
    fn calculate(&self, z: &f64) -> f64 {
        let delta2 = self.delta.powi(2);
        let ds = delta2 + self.sigma2;

        (
            (1.0 + self.delta * z / self.sigma2).powf(
                (1.0 + self.delta * z / self.sigma2) * (- self.sigma2 / ds)
                )
            *
            (1.0 - z / self.delta).powf(
                (z / self.delta - 1.0) * (delta2 / ds)
                )
        ).powi(self.n as i32)
    }
}


impl UpperBound for HoeffdingIneq { 
    fn pr(&self, z: &f64) -> f64 {
        if *z <= self.min_z() {
            1.0
        } else if *z >= self.max_z() {
            0.0
        } else {
            self.calculate(z)
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