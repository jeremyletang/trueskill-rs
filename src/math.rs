
use num::traits::{Float, cast, One, Zero};

fn complementary_error<F: Float + One + Zero>(x: F) -> F {
    let z: F = x.abs();
    let t: F = F::one() / (F::one() + z / cast(2.0).unwrap());
    let r: F = t * (-z * z - cast(1.26551223).unwrap() + t * (cast::<f64, F>(1.00002368).unwrap() +
        t * (cast::<f64, F>(0.37409196).unwrap() + t * (cast::<f64, F>(0.09678418).unwrap() + t * (
            cast::<f64, F>(-0.18628806).unwrap() + t * (cast::<f64, F>(0.27886807).unwrap() + t * (
                cast::<f64, F>(-1.13520398).unwrap() + t * (cast::<f64, F>(1.48851587).unwrap() + t
                * (cast::<f64, F>(-0.82215223).unwrap() + t * cast(0.17087277).unwrap())))))))))
                .exp();
    if x < F::zero() { cast::<f64, F>(2.0).unwrap() - r } else { r }
}

// default to:
// mu: 0, sigma: 1
fn cumulative_distribution<F: Float + One + Zero>(x: F, mu: Option<F>, sigma: Option<F>) -> F {
    let mu_ = if mu.is_some() { mu.unwrap() } else { F::zero() };
    let sigma_ = if sigma.is_some() { sigma.unwrap() } else { F::one() };
    cast::<f64, F>(0.5).unwrap() *
        complementary_error(-(x - mu_) / (sigma_ * cast::<f64, F>(2f64).unwrap().sqrt()))
}

// default to:
// mu: 0, sigma: 1
fn probability_density<F: Float + One + Zero>(x: F, mu: Option<F>, sigma: Option<F>) -> F {
    let mu_ = if mu.is_some() { mu.unwrap() } else { F::zero() };
    let sigma_ = if sigma.is_some() { sigma.unwrap() } else { F::one() };

    (F::one() / (cast::<f64, F>(2f64).unwrap() *
        cast::<f64, F>(::std::f64::consts::PI).unwrap()).sqrt() * sigma_.abs() *
            (-((((x - mu_) / sigma_.abs())).powi(2) / cast::<f64, F>(2f64).unwrap())).exp())
}

fn complementary_error_inv<F: Float + One + Zero>(mut y: F) -> F {
    if y >= cast::<f64, F>(2f64).unwrap() {
        cast::<f64, F>(-100f64).unwrap()
    } else if y <= F::zero() {
        cast::<f64, F>(100f64).unwrap()
    } else {
        let zero_point = y < F::one();
        if !zero_point { y = cast::<f64, F>(2f64).unwrap() - y };
        let t = (cast::<f64, F>(-2f64).unwrap() * y.log(cast::<f64, F>(2f64).unwrap())).sqrt();
        let mut x = cast::<f64, F>(-0.70711).unwrap() * ((cast::<f64, F>(2.30753).unwrap() + t *
            cast::<f64, F>(0.27061).unwrap()) / (F::one() + t * (cast::<f64, F>(0.99229).unwrap() +
            t * cast::<f64, F>(0.04481).unwrap())) - t);

        for _ in 0..2 {
            let err: F = complementary_error(x) - y;
            x = x + err / (cast::<f64, F>(1.12837916709551257).unwrap() * (-(x.powi(2))).exp() - x
                * err);
        }
        if zero_point { x } else { -x }
    }
}

// default to:
// mu: 0, sigma: 1
fn cumulative_distribution_inv<F: Float + One + Zero>(x: F, mu: Option<F>, sigma: Option<F>) -> F {
    let mu_ = if mu.is_some() { mu.unwrap() } else { F::zero() };
    let sigma_ = if sigma.is_some() { sigma.unwrap() } else { F::one() };

    mu_ - sigma_ * cast::<f64, F>(2f64).unwrap().sqrt() *
        complementary_error_inv(cast::<f64, F>(2f64).unwrap() * x)
}

#[test]
fn with_f32() {
    probability_density(10f32, None, None);
    cumulative_distribution(10f32, None, None);
}

#[test]
fn with_f64() {
    probability_density(10f64, None, None);
    cumulative_distribution(10f64, None, None);
}
