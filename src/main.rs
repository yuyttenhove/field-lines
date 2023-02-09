use std::{path::Path, vec};

use glam::DVec2;
use plotlib::{page::Page, repr::Plot, style::LineStyle, view::ContinuousView};

const NEWTON_G: f64 = 1.;
/// Const array of matplotlib colors
const TAB10: [&str; 10] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
];

/// Simple point mass struct
struct PointMass {
    mass: f64,
    pos: DVec2,
}

impl PointMass {
    /// Point mass constructor
    fn new(mass: f64, pos: DVec2) -> Self {
        Self { mass, pos }
    }
}

/// Simple window struct
#[derive(Clone, Copy)]
struct Window {
    anchor: DVec2,
    opposite: DVec2,
}

impl Window {
    fn new(anchor: DVec2, width: DVec2) -> Self {
        Window {
            anchor,
            opposite: anchor + width,
        }
    }

    fn contains(&self, pos: DVec2) -> bool {
        self.anchor.x <= pos.x
            && self.anchor.y <= pos.y
            && pos.x < self.opposite.x
            && pos.y < self.opposite.y
    }
}

/// Calculate the net gravitational acceleration from a slice of point masses
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

/// Compute a field line starting from a given `start_pos`.
/// The Field line is terminated when it exits the viewing `window` or when a stable point is found.
fn field_line(
    start_pos: DVec2,
    point_masses: &[PointMass],
    epsilon: f64,
    window: Window,
) -> Vec<DVec2> {
    let mut line = vec![];
    line.push(start_pos);
    let mut cur_pos = start_pos;
    let mut accel = grav_accel(cur_pos, point_masses);
    let mut accel_norm = accel.length();
    while accel_norm > 1e-6 && window.contains(cur_pos) {
        let dx = -epsilon * accel_norm.min(1.) * accel / accel_norm;
        cur_pos += dx;
        line.push(cur_pos);
        accel = grav_accel(cur_pos, point_masses);
        accel_norm = accel.length();
    }
    line
}

/// Compute the field lines around a given `pos`.
/// For this to work properly, one of the masses from the `point_masses` slice should be located at `pos`.
fn field_lines_around(
    pos: DVec2,
    n: usize,
    point_masses: &[PointMass],
    epsilon: f64,
    window: Window,
) -> Vec<Vec<DVec2>> {
    let mut lines = vec![];
    for i in 0..n {
        let theta = 2. * i as f64 * std::f64::consts::PI / n as f64;
        let start_pos = pos
            + epsilon
                * DVec2 {
                    x: theta.cos(),
                    y: theta.sin(),
                };
        lines.push(field_line(start_pos, point_masses, epsilon, window))
    }
    lines
}

/// Draw the `field_lines' on the `view` with a given `colour`.
fn add_lines(mut view: ContinuousView, field_lines: &[Vec<DVec2>], colour: &str) -> ContinuousView {
    for line in field_lines.iter() {
        let plot = Plot::new(line.iter().map(|v| (v.x, v.y)).collect())
            .line_style(LineStyle::new().colour(colour));
        view = view.add(plot);
    }
    view
}

/// Create a plot the sets of field lines `field_lines` (one set for each mass) with given colorscheme,
/// x/y limits, figure width and savename
fn plot<P: AsRef<Path>>(
    all_field_lines: &[Vec<Vec<DVec2>>],
    colors: &[&str],
    x_range: (f64, f64),
    y_range: (f64, f64),
    width: f64,
    savename: P,
) {
    let mut view = ContinuousView::new();
    for (field_lines, colour) in all_field_lines.iter().zip(colors.iter().cycle()) {
        view = add_lines(view, field_lines, colour);
    }

    view = view
        .x_range(x_range.0, x_range.1)
        .y_range(y_range.0, y_range.1);

    let height = width * (y_range.1 - y_range.0) as f64 / (x_range.1 - x_range.0) as f64;
    Page::single(&view)
        .dimensions(width as u32, height as u32)
        .save(savename)
        .expect("Error saving the plot");
}

fn main() {
    // Construct desired point masses
    let earth = PointMass::new(1., -0.6667 * DVec2::X);
    let moon = PointMass::new(0.02, 0.6667 * DVec2::X);
    let point_masses = vec![earth, moon];

    // generate field lines around each point mass
    let n = 100;
    let epsilon = 0.01;
    let window = Window::new(DVec2 { x: -2., y: -1. }, DVec2 { x: 4., y: 2. });
    let mut all_field_lines = vec![];
    for point_mass in point_masses.iter() {
        all_field_lines.push(field_lines_around(
            point_mass.pos,
            (n as f64 * point_mass.mass.sqrt()) as usize,
            &point_masses,
            epsilon,
            window,
        ));
    }

    plot(
        &all_field_lines,
        &TAB10,
        (window.anchor.x, window.opposite.x),
        (window.anchor.y, window.opposite.y),
        1000.,
        "test.svg",
    );
}
