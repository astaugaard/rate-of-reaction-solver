
use good_lp::{variables, Expression, default_solver, SolverModel, solvers::coin_cbc::CoinCbcSolution, Solution};

// vectors of the form
//
// [reaction_rate, concentration1, .. , concentrationN];

struct Trial {
    rate: f64,
    cons: Vec<f64>
}

type Equations = Vec<Trial>;

#[derive(Debug)]
struct RatesResults {
    k: f64,
    orders: Vec<u64>
}

fn solve(reaction_data: Equations) -> Option<RatesResults> {
    let eqs = reaction_data.len();
    let num_orders = reaction_data.get(0)?.cons.len();

    variables! {
        vars:
            errorsp[eqs] >= 0;
            errorsm[eqs] >= 0;
            orders[num_orders] (integer);
            pk >= 0;
        };

    let expr:Expression =
        errorsp.iter().map(|expr| {
        let expr:Expression = (*expr).into();
        expr
        }).sum::<Expression>() +
        errorsm.iter().map(|expr| {
            let expr:Expression = (*expr).into();
            expr}).sum::<Expression>();

    let mut solution = vars.minimise(expr).using(default_solver);

    for (trial_num,trial) in reaction_data.iter().enumerate() {
        let error:Expression = Expression::from(errorsp[trial_num]) - Expression::from(errorsm[trial_num]);
        let rate:Expression = trial.rate.ln().into();
        let pk:Expression = pk.into();
        let con = trial.cons.iter()
            .enumerate()
            .map(|(o,concentration)| -> Expression {
                let order:Expression = orders[o].into();
                let lncon= concentration.ln();
                order * lncon})
            .fold(error + pk, |a:Expression,b| {
                let expr = a;
                expr + b}).eq(rate);
        solution = solution.with(con);
    }

    let result:CoinCbcSolution = solution.solve().unwrap();

    let orders:Vec<u64> = (0..num_orders).map(|i| result.value(orders[i]) as u64).collect();
    let k = result.value(pk).exp();

    Some(RatesResults { k: k, orders: orders })
}

fn main() {
    let t1 = Trial {rate: 0.136, cons: vec![0.420,0.122]};
    let t2 = Trial {rate: 0.0339, cons: vec![0.210,0.122]};
    let t3 = Trial {rate: 0.0678, cons: vec![0.21,0.244]};
    let t4 = Trial {rate: 0.0339, cons: vec![0.105,0.488]};

    println!("{:#?}", solve(vec![t1,t2,t3,t4]));

}
