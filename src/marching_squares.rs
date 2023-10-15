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

pub(crate) fn pick_bits(f: &[f32], shape: Shape, pos: (isize, isize), threshold: f32) -> u8 {
    (f[shape.idx(pos.0, pos.1)] > threshold) as u8
        | (((f[shape.idx(pos.0 + 1, pos.1)] > threshold) as u8) << 1)
        | (((f[shape.idx(pos.0 + 1, pos.1 + 1)] > threshold) as u8) << 2)
        | (((f[shape.idx(pos.0, pos.1 + 1)] > threshold) as u8) << 3)
}

/// LINE_WIDTH won't work well with cargo fmt
const LW: f32 = 0.;

pub(crate) const CELL_BORDER_BUFFER: [[f32; 4]; 6] = [
    [1., 0., -1., 0.],
    [0., 1., 0., -1.],
    [-1., 0., 0., -1.],
    [0., -1., 1., 0.],
    [1., 0., 0., 1.],
    [0., 1., -1., 0.],
];

/// Index into CELL_BORDER_BUFFER
pub(crate) fn cell_border_index(bits: u8) -> Option<&'static [f32]> {
    let idx = match bits {
        1 | 14 => 2,
        2 | 13 => 3,
        4 | 11 => 4,
        8 | 7 => 5,
        3 | 12 => 0,
        9 | 6 => 1,
        5 => return Some(&[-1., 0., 0., -1., 1., 0., 0., 1.]),
        10 => return Some(&[0., -1., 1., 0., 0., 1., -1., 0.]),
        _ => return None,
    };
    Some(&CELL_BORDER_BUFFER[idx])
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
