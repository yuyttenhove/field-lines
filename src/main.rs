use std::{path::Path, vec};

use glam::DVec2;
use plotlib::{page::Page, repr::Plot, style::LineStyle, view::ContinuousView};

const NEWTON_G: f64 = 1.;
// pyplot colors
const TAB10: [&str; 10] = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];

struct PointMass {
    mass: f64,
    pos: DVec2,
}

impl PointMass {
    fn new(mass: f64, pos: DVec2) -> Self {
        Self { mass, pos }
    }
}

fn grav_accel(pos: DVec2, point_masses: &[PointMass]) -> DVec2 {
    fn grav_accel_single(pos: DVec2, pm: &PointMass) -> DVec2 {
        let dx = pm.pos - pos;
        NEWTON_G * pm.mass / dx.length_squared() * dx.normalize()
    }
    point_masses
        .iter()
        .map(|pm| grav_accel_single(pos, pm))
        .sum()
}

fn field_line(
    start_pos: DVec2,
    point_masses: &[PointMass],
    alpha: f64,
    max_dist: f64,
) -> Vec<DVec2> {
    let mut line = vec![];
    line.push(start_pos);
    let mut cur_pos = start_pos;
    let mut accel = grav_accel(cur_pos, point_masses);
    let mut accel_norm = accel.length();
    while accel_norm > 1e-8
        && (cur_pos - start_pos).length_squared() < max_dist * max_dist
    {
        let dx = -alpha * accel_norm.min(1.) * accel / accel_norm;
        cur_pos += dx;
        line.push(cur_pos);
        accel = grav_accel(cur_pos, point_masses);
        accel_norm = accel.length();
    }
    line
}

fn field_lines_around(pos: DVec2, n: usize, point_masses: &[PointMass], alpha: f64, max_dist: f64) -> Vec<Vec<DVec2>> {
    let mut lines = vec![];
    for i in 0..n {
        let theta = 2. * i as f64 * std::f64::consts::PI / n as f64;
        let start_pos = pos + alpha * DVec2 { x: theta.cos(), y: theta.sin() };
        lines.push(field_line(start_pos, point_masses, alpha, max_dist))
    }
    lines
}

fn add_lines(v: ContinuousView, field_lines: &[Vec<DVec2>], colour: &str) -> ContinuousView {
    let mut v = v;
    for line in field_lines.iter() {
        let plot =
            Plot::new(line.iter().map(|v| (v.x, v.y)).collect()).line_style(LineStyle::new().colour(colour));
        v = v.add(plot);
    }
    v
}

fn plot<P: AsRef<Path>>(
    field_lines: &[Vec<Vec<DVec2>>],
    colors: &[&str],
    x_range: (f64, f64),
    y_range: (f64, f64),
    width: f64,
    savename: P,
) {
    let mut v = ContinuousView::new();
    for (lines, colour) in field_lines.iter().zip(colors.iter().cycle()) {
        v = add_lines(v, lines, colour);
    }

    v = v
        .x_range(x_range.0, x_range.1)
        .y_range(y_range.0, y_range.1);

    let height = width * (y_range.1 - y_range.0) as f64 / (x_range.1 - x_range.0) as f64;
    Page::single(&v)
        .dimensions(width as u32, height as u32)
        .save(savename)
        .expect("Error saving the plot");
}

fn main() {
    let point_masses = vec![PointMass::new(1., -DVec2::X), PointMass::new(0.1, DVec2::X)];

    // generate field lines around each point mass
    let n = 90;
    let alpha = 0.01;
    let max_dist = 6.;
    let mut field_lines = vec![];
    for pm in point_masses.iter() {
        field_lines.push(field_lines_around(pm.pos, (n as f64 * pm.mass.sqrt()) as usize, &point_masses, alpha, max_dist));
    }

    plot(&field_lines,  &TAB10, (-3., 3.), (-3., 3.), 1000., "test.svg");
}
