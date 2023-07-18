
use pyo3::prelude::*;
use eqsolver::ODESolver;

#[pymodule]
fn rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(is_prime, m)?)?;
    m.add_function(wrap_pyfunction!(is_prime2, m)?)?;
    m.add_function(wrap_pyfunction!(ode_example, m)?)?;
    Ok(())
}

  
#[pyfunction] 
fn is_prime(num: u32) -> bool {
    match num {
        0 | 1 => false,
        _ => {
            let limit = (num as f32).sqrt() as u32; 

            (2..=limit).any(|i| num % i == 0) == false
        }
    }
}
#[pyfunction] 
fn is_prime2(num: u32) -> bool {
    match num {
        0 | 1 => false,
        _ => {
            let limit = (num as f32).sqrt() as u32; 

            (2..=limit).any(|i| num % i == 0) == false
        }
    }
}
#[pyfunction]
fn ode_example(x0: f64, y0: f64,x_end: f64, step_size: f64) -> f64{
    let f = |t: f64, y: f64| t /y/y;
    let solver = ODESolver::new(f, x0,y0, step_size);
    let solution: f64  = solver.solve(x_end).unwrap();
    return solution;
}