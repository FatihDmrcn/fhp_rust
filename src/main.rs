use std::collections::HashMap;
use std::f32::consts::PI;
use std::fmt;
use std::hash::Hash;
use std::ops::Index;
use rayon::prelude::*;
use rand::random;
use rand::seq::SliceRandom;

use plotters::prelude::*;
use indicatif::ProgressBar;



#[derive(Debug, PartialEq, Eq, Hash, Clone)]
enum CellType { Air, Wall }

impl Default for CellType { fn default() -> Self { CellType::Air } }

#[repr(u64)]
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
enum Dir {
    R   = 0,
    TR  = 1,
    TL  = 2,
    L   = 3,
    BL  = 4,
    BR  = 5,
}

struct Scatter {}

impl Scatter {
    fn v2a() -> [u8; 6] { vec![[0,1,0,0,1,0], [1,0,0,1,0,0]].choose(&mut rand::thread_rng()).unwrap().clone() }
    fn v2b() -> [u8; 6] { vec![[0,0,1,0,0,1], [1,0,0,1,0,0]].choose(&mut rand::thread_rng()).unwrap().clone() }
    fn v2c() -> [u8; 6] { vec![[0,0,1,0,0,1], [0,1,0,0,1,0]].choose(&mut rand::thread_rng()).unwrap().clone() }
    fn v3a() -> [u8; 6] { [1,0,1,0,1,0] }
    fn v3b() -> [u8; 6] { [0,1,0,1,0,1] }
    fn v4a() -> [u8; 6] { vec![[1,0,1,1,0,1], [1,1,0,1,1,0]].choose(&mut rand::thread_rng()).unwrap().clone() }
    fn v4b() -> [u8; 6] { vec![[0,1,1,0,1,1], [1,1,0,1,1,0]].choose(&mut rand::thread_rng()).unwrap().clone() }
    fn v4c() -> [u8; 6] { vec![[0,1,1,0,1,1], [1,0,1,1,0,1]].choose(&mut rand::thread_rng()).unwrap().clone() }
}


impl Dir {
    const VALUES: [Self; 6] = [Self::R,Self::TR,Self::TL,Self::L,Self::BL,Self::BR];
    fn opposite(&self) -> Self {
        match *self {
            Dir::R  => Dir::L,
            Dir::TR => Dir::BL,
            Dir::TL => Dir::BR,
            Dir::L  => Dir::R,
            Dir::BL => Dir::TR,
            Dir::BR => Dir::TL,
        }
    }
}

struct Universe{
    rows: u64,
    cols: u64,
    cells: Vec<Cell>,
}

impl Universe {
    fn new(rows: &u64, cols: &u64) -> Self{
        Self {
            rows: rows.clone(),
            cols: cols.clone(),
            cells: Self::__init_cells(rows, cols),
        }
    }
    fn __init_cells(rows: &u64, cols: &u64) -> Vec<Cell> {
        /*
        #EXEMPLARY GRID#
        00  01  02  03
          04  05  06  07  
        08  09  10  11
          12  13  14  15
        */
        let mut cells: Vec<Cell> = vec![];
        for _r in 0..*rows {
            for _c in 0..*cols {
                let _idx = _r*cols + _c;
                let _last_col = cols-1;
                let _last_row = rows-1;
                let _mod_row = _r%2;
                let mut adjacency = HashMap::new();
                if _c != _last_col {adjacency.insert(Dir::R, _idx+1);}
                if _c != 0 {adjacency.insert(Dir::L, _idx-1);}
                if _r != _last_row {
                    if _mod_row==0 {
                        adjacency.insert(Dir::BR, _idx+cols);
                        if _c != 0 {adjacency.insert(Dir::BL, _idx+cols-1);}
                    }
                    else {
                        adjacency.insert(Dir::BL, _idx+cols);
                        if _c != _last_col {adjacency.insert(Dir::BR, _idx+cols+1);}
                    }
                }
                if _r != 0 {
                    if _mod_row==0 {
                        adjacency.insert(Dir::TR, _idx-cols);
                        if _c != 0 {adjacency.insert(Dir::TL, _idx-(cols+1));}
                    }
                    else {
                        adjacency.insert(Dir::TL, _idx-cols);
                        if _c != _last_col {adjacency.insert(Dir::TR, _idx-(cols-1));}
                    }
                }
                let mut cell_type = CellType::Wall;
                let mut neumann_bc: Option<[u8;6]> = None;
                //
                let walls = _r == 0 || _r == _last_row;
                if !walls {
                    cell_type = CellType::Air;
                    if _c == 0 {neumann_bc = Some([1,1,0,0,0,1])}
                }
                cells.push( Cell::new(_idx.to_owned(), cell_type, Some(adjacency), neumann_bc, Some(_r), Some(_c)) );
            }
        }
        cells
    }
    fn next_step(_cells: &Vec<Cell>) -> Vec<Cell> {
        let mut new_cells: Vec<Cell> = vec![];
        _cells.into_par_iter().map(|cell| {
            if cell.cell_type != CellType::Wall {
            // 0th - IF APPLICABLE - NEUMANN BC, OR ELSE
            let mut _vel: [u8; 6] = cell.nbc.unwrap_or_else(|| {
                let mut composed_vel: [u8; 6] = [0,0,0,0,0,0];
                // 1st - FLUX
                for _dir in Dir::VALUES {
                    let _idx_dir = _dir.clone() as usize;
                    let _idx_dir_opposed = _dir.opposite() as usize;
                    let _vel_dir = cell.adj.get(&_dir).map(|_neighbor_id| {
                        let mut _vel: u8;
                        let _neighbor_at_dir = _cells.index(*_neighbor_id as usize); 
                        if _neighbor_at_dir.cell_type == CellType::Wall { _vel = cell.vel[_idx_dir]; }
                        else { _vel = _neighbor_at_dir.vel[_idx_dir_opposed]; }
                        _vel
                    }
                    );
                    composed_vel[_idx_dir_opposed] = _vel_dir.unwrap_or(0);
                }
                // 2nd - APPLY SCATTERING
                match composed_vel {
                    [0,0,1,0,0,1] => { composed_vel = Scatter::v2a(); },
                    [0,1,0,0,1,0] => { composed_vel = Scatter::v2b(); },
                    [1,0,0,1,0,0] => { composed_vel = Scatter::v2c(); },
                    [0,1,0,1,0,1] => { composed_vel = Scatter::v3a(); },
                    [1,0,1,0,1,0] => { composed_vel = Scatter::v3b(); },
                    [0,1,1,0,1,1] => { composed_vel = Scatter::v4a(); },
                    [1,0,1,1,0,1] => { composed_vel = Scatter::v4b(); },
                    [1,1,0,1,1,0] => { composed_vel = Scatter::v4c(); },
                    _ => {}
                }
                composed_vel
                });
            cell.copy(_vel)}
            else {cell.clone()}
        }).collect_into_vec(&mut new_cells);
        // 3rd - ENSURE ASCENDING SORTING
        new_cells.sort_by(|a, b| a.id.cmp(&b.id));
        new_cells
    }
}

#[derive(Default, Debug, PartialEq, Clone)]
struct Cell{
    id: u64,
    cell_type: CellType,
    adj: HashMap<Dir,u64>,
    nbc: Option<[u8; 6]>,
    vel: [u8; 6],
    row: u64,
    col: u64,
}

impl Cell {
    fn new(
        id: u64, 
        cell_type: CellType, 
        adjacency: Option<HashMap<Dir,u64>>, 
        neumann_bc: Option<[u8; 6]>, 
        row: Option<u64>, 
        col: Option<u64>) -> Cell {

        let mut _vel: [u8; 6] = [0,0,0,0,0,0];
        // if cell_type == CellType::Air { _vel = neumann_bc.unwrap_or([0,0,0,0,0,0]) }
        if cell_type == CellType::Air {
            _vel = neumann_bc.unwrap_or_else(|| {
                let mut _vel_rand: [u8; 6] = _vel.clone();
                for _dir in Dir::VALUES {_vel_rand[_dir as usize] = random::<bool>() as u8;}
                _vel_rand
            })
        }
        Cell {
            id,
            cell_type,
            adj: adjacency.unwrap_or(HashMap::new()),
            nbc: neumann_bc,
            vel: _vel,
            row: row.unwrap_or(0),
            col: col.unwrap_or(0),
        }
    }
    fn copy(&self, new_vel: [u8; 6]) -> Cell {
        Cell {
            id: self.id.clone(),
            cell_type: self.cell_type.clone(),
            adj: self.adj.clone(),
            nbc: self.nbc.clone(),
            vel: new_vel,
            row: self.row.clone(),
            col: self.col.clone(),
        }
    }
    fn total_vel(&self) -> u8 {
        self.vel.into_iter().sum()
    }
    fn vel_h(&self) -> f64 {
        let v = self.vel.map(f64::from);
        let a = (PI/3.).cos() as f64;
        v[0] + a*v[1] - a*v[2] - v[3] - a*v[4] + a*v[5]
    }
    fn vel_v(&self) -> f64 {
        let v = self.vel.map(f64::from);
        let a = (PI/3.).sin() as f64;
        a*v[1] + a*v[2] - a*v[4] - a*v[5]
    } 
}

impl fmt::Display for Universe {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _print: Vec<String> = vec![];
        let mut _row_str: Vec<String> = vec![];
        for _r in 0..self.rows {
            let mut _col_str: Vec<String> = vec![];
            for _c in 0..self.cols {
                let _idx = _r*self.cols+_c;
                _col_str.push(format!("{:<4}", self.cells.get(_idx as usize).unwrap().total_vel()))
            }
            let _row: String = _col_str.join("");
            let _spacer: &str;
            if _r%2==0 {_spacer="";}
            else {_spacer="  ";}
            _row_str.push(format!("{}{}", _spacer, _row));
        }
        write!(f, "{}", _row_str.join("\n"))
    }
}

impl fmt::Display for Cell {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Cell ID {}\tCell Type{}\tParticles {:?}", self.id, self.cell_type, self.vel)
    }
}

impl fmt::Display for CellType {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            CellType::Air => write!(f, "Air"),
            CellType::Wall => write!(f, "Wall"),
        }
    }
}

impl fmt::Display for Dir {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Dir::L => write!(f, " L"),
            Dir::R => write!(f, " R"),
            Dir::TL => write!(f, "TL"),
            Dir::TR => write!(f, "TR"),
            Dir::BL => write!(f, "BL"),
            Dir::BR => write!(f, "BR"),
        }
    }
}


fn main() {
    let rows = 40;
    let cols = 80; //84;
    let mut universe = Universe::new(&rows, &cols);
    //
    let r = 4;
    let d = 2*r;
    let w = cols*d+r;
    let h = d*rows-(r/2);
    let t_total = 201;
    let area = BitMapBackend::gif("misc/animated.gif", (w as u32, h as u32), 10/*ms*/).unwrap().into_drawing_area();
    let pb = ProgressBar::new(t_total);
    for _ in 0..=t_total {
        area.fill(&WHITE).unwrap();
        for cell in universe.cells.iter() {
            let _r: i32 = cell.row as i32;
            let _c: i32 = cell.col as i32;
            let mut _coords: (i32, i32) = (r as i32+_c*d as i32, (r as i32/2)+_r*d as i32);
            if _r%2!=0 { _coords = (d as i32+_c*d as i32, (r as i32/2)+_r*d as i32); }
            let mut _color = ViridisRGBA::get_color_normalized(cell.vel_h(),-2.,2.);
            if cell.cell_type == CellType::Wall {_color = BLACK.into();}
            let _style = ShapeStyle {color: _color, filled: true, stroke_width: 0};
            area.draw(&Circle::new(_coords, r as i32, _style)).unwrap();
        }
        area.present().unwrap();
        pb.inc(1);
        universe.cells = Universe::next_step(&universe.cells);
        }
    pb.finish_with_message("done");
}