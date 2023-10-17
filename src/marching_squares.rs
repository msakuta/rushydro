//! Marching squares algorithm
//!
//! In marching squares, the coordinates specification is very important.
//! Each vertex of a cell is indexed like below.
//! Note that the coordinates are in cartesian, meaning positive Y is pointing up,
//! which is the opposite of most screen coordinates.
//!
//!       ^
//!  3 ---+--- 2
//!  |    |    |
//!  +----+----+--->
//!  |    |    |
//!  0----+----1

pub(crate) trait Idx {
    fn idx(&self, x: isize, y: isize) -> usize;
}

pub(crate) type Shape = (isize, isize);

impl Idx for Shape {
    fn idx(&self, x: isize, y: isize) -> usize {
        let (width, height) = self;
        ((x + width) % width + (y + height) % height * width) as usize
    }
}

/// Returns whether each vertex of the cell is above threshold in a packed byte.
/// LSB is the 0th vertex.
pub(crate) fn pick_bits(f: &[f32], shape: Shape, pos: (isize, isize), threshold: f32) -> u8 {
    (f[shape.idx(pos.0, pos.1)] > threshold) as u8
        | (((f[shape.idx(pos.0 + 1, pos.1)] > threshold) as u8) << 1)
        | (((f[shape.idx(pos.0 + 1, pos.1 + 1)] > threshold) as u8) << 2)
        | (((f[shape.idx(pos.0, pos.1 + 1)] > threshold) as u8) << 3)
}

pub(crate) fn pick_values(f: &[f32], shape: Shape, pos: (isize, isize)) -> [f32; 4] {
    [
        f[shape.idx(pos.0, pos.1)],
        f[shape.idx(pos.0 + 1, pos.1)],
        f[shape.idx(pos.0 + 1, pos.1 + 1)],
        f[shape.idx(pos.0, pos.1 + 1)],
    ]
}

/// LINE_WIDTH won't work well with cargo fmt
const LW: f32 = 0.;

#[allow(dead_code)]
pub(crate) fn cell_border_index(bits: u8) -> Option<&'static [f32]> {
    Some(match bits {
        1 | 14 => &[-1., 0., 0., -1.],
        2 | 13 => &[0., -1., 1., 0.],
        4 | 11 => &[1., 0., 0., 1.],
        8 | 7 => &[0., 1., -1., 0.],
        3 | 12 => &[1., 0., -1., 0.],
        9 | 6 => &[0., 1., 0., -1.],
        5 => &[-1., 0., 0., -1., 1., 0., 0., 1.],
        10 => &[0., -1., 1., 0., 0., 1., -1., 0.],
        _ => return None,
    })
}

pub(crate) fn cell_border_interpolated(bits: u8, values: [f32; 4]) -> Option<([f32; 8], usize)> {
    let factor = |f0: f32, f1| (1. - (f0 + f1)) / (-f0 + f1);
    let arr = match bits {
        1 | 14 => {
            let x = factor(values[0], values[1]);
            let y = factor(values[0], values[3]);
            [-1., y, x, -1.]
        }
        2 | 13 => {
            let x = factor(values[0], values[1]);
            let y = factor(values[1], values[2]);
            [x, -1., 1., y]
        }
        4 | 11 => {
            let x = factor(values[3], values[2]);
            let y = factor(values[1], values[2]);
            [1., y, x, 1.]
        }
        8 | 7 => {
            let x = factor(values[3], values[2]);
            let y = factor(values[0], values[3]);
            [x, 1., -1., y]
        }
        3 | 12 => {
            let y0 = factor(values[0], values[3]);
            let y1 = factor(values[1], values[2]);
            [1., y1, -1., y0]
        }
        9 | 6 => {
            let x0 = factor(values[0], values[1]);
            let x1 = factor(values[3], values[2]);
            [x1, 1., x0, -1.]
        }
        5 => {
            let x0 = factor(values[0], values[1]);
            let y0 = factor(values[0], values[3]);
            let x1 = factor(values[3], values[2]);
            let y1 = factor(values[1], values[2]);
            return Some(([-1., y0, x0, -1., 1., y1, x1, 1.], 8));
        }
        10 => {
            let x0 = factor(values[0], values[1]);
            let y0 = factor(values[1], values[2]);
            let x1 = factor(values[3], values[2]);
            let y1 = factor(values[0], values[3]);
            return Some(([x0, -1., 1., y0, x1, 1., -1., y1], 8));
        }
        _ => return None,
    };
    let mut ret = [0.; 8];
    ret[0..4].copy_from_slice(&arr);
    Some((ret, 4))
}

pub(crate) fn cell_polygon_interpolated(bits: u8, values: [f32; 4]) -> Option<([f32; 10], usize)> {
    let factor = |f0: f32, f1| (1. - (f0 + f1)) / (-f0 + f1);
    let f = |i, j| factor(values[i], values[j]);
    let mut ret = [0.; 10];
    let num = match bits {
        1 => {
            ret[..6].copy_from_slice(&[-1., f(0, 3), f(0, 1), -1., -1., -1.]);
            6
        }
        14 => {
            ret[..10].copy_from_slice(&[-1., f(0, 3), f(0, 1), -1., 1., -1., 1., 1., -1., 1.]);
            10
        }
        2 => {
            ret[..6].copy_from_slice(&[f(0, 1), -1., 1., f(1, 2), 1., -1.]);
            6
        }
        13 => {
            ret[..10].copy_from_slice(&[f(0, 1), -1., 1., f(1, 2), 1., 1., -1., 1., -1., -1.]);
            10
        }
        4 => {
            ret[..6].copy_from_slice(&[1., f(1, 2), f(3, 2), 1., 1., 1.]);
            6
        }
        11 => {
            ret[..10].copy_from_slice(&[1., f(1, 2), f(3, 2), 1., -1., 1., -1., -1., 1., -1.]);
            10
        }
        8 => {
            ret[..6].copy_from_slice(&[f(3, 2), 1., -1., f(0, 3), -1., 1.]);
            6
        }
        7 => {
            ret[..10].copy_from_slice(&[f(3, 2), 1., -1., f(0, 3), -1., -1., 1., -1., 1., 1.]);
            10
        }
        3 => {
            ret[..8].copy_from_slice(&[1., f(1, 2), -1., f(0, 3), -1., -1., 1., -1.]);
            8
        }
        12 => {
            ret[..8].copy_from_slice(&[1., f(1, 2), -1., f(0, 3), -1., 1., 1., 1.]);
            8
        }
        9 => {
            ret[..8].copy_from_slice(&[f(3, 2), 1., f(0, 1), -1., -1., -1., -1., 1.]);
            8
        }
        6 => {
            ret[..8].copy_from_slice(&[f(3, 2), 1., f(0, 1), -1., 1., -1., 1., 1.]);
            8
        }
        5 => {
            let x0 = factor(values[0], values[1]);
            let y0 = factor(values[0], values[3]);
            let x1 = factor(values[3], values[2]);
            let y1 = factor(values[1], values[2]);
            ret[..8].copy_from_slice(&[-1., y0, x0, -1., 1., y1, x1, 1.]);
            8
        }
        10 => {
            let x0 = factor(values[0], values[1]);
            let y0 = factor(values[1], values[2]);
            let x1 = factor(values[3], values[2]);
            let y1 = factor(values[0], values[3]);
            ret[..8].copy_from_slice(&[x0, -1., 1., y0, x1, 1., -1., y1]);
            8
        }
        _ => return None,
    };
    Some((ret, num))
}

/// buffer for vertex shader, use with SliceFlatExt::flat()
#[allow(dead_code)]
pub(crate) const CELL_POLYGON_BUFFER: [[f32; 8]; 7] = [
    [1., 1., -1., 1., -1., -1., 1., -1.],
    [1., LW, -1., LW, -1., -LW, 1., -LW],
    [LW, 1., -LW, 1.0, -LW, -1., LW, -1.],
    [-1., -LW, -LW, -1., LW, -1., -1., LW],
    [LW, -1., 1., -LW, 1., LW, -LW, -1.],
    [1., LW, LW, 1., -LW, 1., 1., -LW],
    [-LW, 1., -1., LW, -1., -LW, LW, 1.],
];

/// Index into CELL_POLYGON_BUFFER
#[allow(dead_code)]
pub(crate) fn cell_polygon_index(bits: u8) -> u8 {
    match bits {
        1 | 14 => 12,
        2 | 13 => 16,
        4 | 11 => 20,
        8 | 7 => 24,
        3 | 12 => 4,
        9 | 6 => 8,
        _ => 0,
    }
}

/// Whether the pixel is a border
pub(crate) fn border_pixel(idx: u8) -> bool {
    match idx {
        0 => false,
        1..=14 => true,
        15 => false,
        _ => panic!("index must be in 0-15!"),
    }
}

#[test]
fn test_bits() {
    assert_eq!(pick_bits(&[0., 0., 0., 0.], (2, 2), (0, 0), 0.5), 0);
    assert_eq!(pick_bits(&[1., 0., 0., 0.], (2, 2), (0, 0), 0.5), 1);
    assert_eq!(pick_bits(&[0., 1., 0., 0.], (2, 2), (0, 0), 0.5), 2);
    assert_eq!(pick_bits(&[0., 0., 1., 0.], (2, 2), (0, 0), 0.5), 8);
    assert_eq!(pick_bits(&[0., 0., 0., 1.], (2, 2), (0, 0), 0.5), 4);
}
