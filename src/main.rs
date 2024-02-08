use hoeffding_ineq::polynominal::PolynominalIneq;

fn main() {
    println!("確率不等式計算プログラム");

    // 個々の値をセットする場合のコード例
    // let sigma2 = 4.0;
    // let delta = 8.0;
    // let n = 3;
    // let pie = PolynominalIneq::new(n, sigma2, delta).unwrap();
    // println!("{:?}", pie);
    
    // let path_test = format!("test/test_polynominal_{n}.xlsx");
    // pie.compare_hoeffding_to_excel(&path_test).unwrap();
    // let path_overview = format!("test/overview_polynominal_{n}.xlsx");
    // pie.overview_with_root_to_excel(&path_overview).unwrap();

    data_for_jima().unwrap();
}

// 論文投稿データ作成用
fn data_for_jima() -> Result<(), Box<dyn std::error::Error>> {
    // 理論値との比較
    let sigma2 = 4.0;
    let delta = 5.0;
    let ns = [2, 5, 10];
    for n in ns {
        let hoef_mitsu = PolynominalIneq::new(n, sigma2, delta)?;
        let path_compare = format!("test/compare_hoeffding_n={}.xlsx", n);
        let path_overview = format!("test/overview_n={}.xlsx", n);
        hoef_mitsu.compare_hoeffding_to_excel(&path_compare)?;
        hoef_mitsu.overview_with_root_to_excel(&path_overview)?;
        let z_sig2 = (sigma2 * (n as f64) / (n as f64).powi(2)).sqrt() * 2.0;
        let z_round = (z_sig2 * 10.0).round() / 10.0;
        let path_root_sig2 = format!("test/root_n={}_z=sig2.xlsx", n);
        let path_root_round = format!("test/root_n={}_10z={}.xlsx", n, z_round * 10.0);
        hoef_mitsu.info_each_root_to_excel(z_sig2, &path_root_sig2)?;
        hoef_mitsu.info_each_root_to_excel(z_round, &path_root_round)?;
    }

    println!("Calculated");
    Ok(())
}
