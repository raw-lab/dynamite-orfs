use std::ffi::OsStr;
use std::fs::File;
use std::io::ErrorKind;
use std::io::{self, BufRead, Lines};
use std::path::PathBuf;
use std::str::FromStr;

extern crate thiserror;
use thiserror::Error;

extern crate strum;
use strum::EnumCount;
use strum_macros::{EnumCount, EnumIter};

use crate::fgs::dna::{ACGT, BI_ACGT, CG_MAX, CG_MIN, TRI_ACGT};

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, EnumCount, EnumIter)]
pub enum State {
    S = 0,
    E = 1,
    R = 2,
    Sr = 3,
    Er = 4,
    M1 = 5,
    M2 = 6,
    M3 = 7,
    M4 = 8,
    M5 = 9,
    M6 = 10,
    M1r = 11,
    M2r = 12,
    M3r = 13,
    M4r = 14,
    M5r = 15,
    M6r = 16,
    I1 = 17,
    I2 = 18,
    I3 = 19,
    I4 = 20,
    I5 = 21,
    I6 = 22,
    I1r = 23,
    I2r = 24,
    I3r = 25,
    I4r = 26,
    I5r = 27,
    I6r = 28,
}

impl<T> std::ops::Index<State> for [T] {
    type Output = T;

    fn index(&self, idx: State) -> &T {
        &self[idx as usize]
    }
}

impl<T> std::ops::IndexMut<State> for [T] {
    fn index_mut(&mut self, idx: State) -> &mut T {
        &mut self[idx as usize]
    }
}

const NUM_TRANSITIONS: usize = 14;

#[derive(Default)]
pub struct Transition {
    pub mm: f64,
    pub mi: f64,
    pub md: f64,
    pub ii: f64,
    pub im: f64,
    pub dd: f64,
    pub dm: f64,
    pub ge: f64,
    pub gg: f64,
    pub er: f64,
    pub rs: f64,
    pub rr: f64,
    pub es: f64,
    pub es1: f64,
}

pub const PERIOD: usize = 6;
const WINDOW: usize = 61;

#[derive(Default)]
pub struct Global {
    pub pi: [f64; State::COUNT],
    pub tr: Transition,
    pub tr_ii: [[f64; ACGT]; ACGT],
    pub tr_mi: [[f64; ACGT]; ACGT],
}

pub struct Local {
    pub e_m: [[[f64; ACGT]; BI_ACGT]; PERIOD],
    pub e_m1: [[[f64; ACGT]; BI_ACGT]; PERIOD],

    pub tr_rr: [[f64; ACGT]; ACGT],

    pub tr_s: [[f64; TRI_ACGT]; WINDOW],
    pub tr_e: [[f64; TRI_ACGT]; WINDOW],
    pub tr_s1: [[f64; TRI_ACGT]; WINDOW],
    pub tr_e1: [[f64; TRI_ACGT]; WINDOW],

    pub dist_s: [f64; PERIOD],
    pub dist_e: [f64; PERIOD],
    pub dist_s1: [f64; PERIOD],
    pub dist_e1: [f64; PERIOD],
}

#[derive(Error, Debug)]
pub enum TrainingDataError {
    #[error("incomplete training file")]
    IncompleteTrainingFile,
    #[error("unknown transition state")]
    UnknownTransitionState,
    #[error("could not read training file {0}")]
    Io(PathBuf, #[source] io::Error),
    #[error("malformed value '{2}' in section {1} of {0}")]
    MalformedValue(PathBuf, String, String, #[source] std::num::ParseFloatError),
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(Box<Global>, Vec<Local>), TrainingDataError> {
    let mut global: Box<Global> = Box::new(Default::default());
    let mut locals: Vec<Local> = (0..CG_MAX - CG_MIN)
        .map(|_| Local {
            e_m: [[[0.0; ACGT]; BI_ACGT]; PERIOD],
            e_m1: [[[0.0; ACGT]; BI_ACGT]; PERIOD],
            tr_rr: [[0.0; ACGT]; ACGT],
            tr_s: [[0.0; TRI_ACGT]; WINDOW],
            tr_e: [[0.0; TRI_ACGT]; WINDOW],
            tr_s1: [[0.0; TRI_ACGT]; WINDOW],
            tr_e1: [[0.0; TRI_ACGT]; WINDOW],

            dist_s: [0.0; PERIOD],
            dist_e: [0.0; PERIOD],
            dist_s1: [0.0; PERIOD],
            dist_e1: [0.0; PERIOD],
        })
        .collect();

    read_transitions(&mut global, train_dir.join(filename))?;
    read_m_transitions(&mut locals, train_dir.join("gene"))?;
    read_m1_transitions(&mut locals, train_dir.join("rgene"))?;
    read_noncoding(&mut locals, train_dir.join("noncoding"))?;
    read_start(&mut locals, train_dir.join("start"))?;
    read_stop(&mut locals, train_dir.join("stop"))?;
    read_start1(&mut locals, train_dir.join("stop1"))?; // keep FGS naming scheme
    read_stop1(&mut locals, train_dir.join("start1"))?; // keep FGS naming scheme
    read_pwm(&mut locals, train_dir.join("pwm"))?;

    Ok((global, locals))
}

fn lines_from_bytes(bytes: &'static [u8]) -> Lines<Box<dyn BufRead>> {
    let b: Box<dyn BufRead> = Box::new(bytes);
    b.lines()
}

fn lines_from_file(filename: &PathBuf) -> Result<Lines<Box<dyn BufRead>>, TrainingDataError> {
    match File::open(filename) {
        Ok(file) => {
            let b: Box<dyn BufRead> = Box::new(io::BufReader::new(file));
            Ok(b.lines())
        }
        Err(e) if e.kind() == ErrorKind::NotFound => {
            match filename.file_name().and_then(OsStr::to_str) {
                Some("454_10") => Ok(lines_from_bytes(include_bytes!("train/454_10"))),
                Some("454_30") => Ok(lines_from_bytes(include_bytes!("train/454_30"))),
                Some("454_5") => Ok(lines_from_bytes(include_bytes!("train/454_5"))),
                Some("complete") => Ok(lines_from_bytes(include_bytes!("train/complete"))),
                Some("gene") => Ok(lines_from_bytes(include_bytes!("train/gene"))),
                Some("illumina_1") => Ok(lines_from_bytes(include_bytes!("train/illumina_1"))),
                Some("illumina_10") => Ok(lines_from_bytes(include_bytes!("train/illumina_10"))),
                Some("illumina_5") => Ok(lines_from_bytes(include_bytes!("train/illumina_5"))),
                Some("noncoding") => Ok(lines_from_bytes(include_bytes!("train/noncoding"))),
                Some("pwm") => Ok(lines_from_bytes(include_bytes!("train/pwm"))),
                Some("rgene") => Ok(lines_from_bytes(include_bytes!("train/rgene"))),
                Some("sanger_10") => Ok(lines_from_bytes(include_bytes!("train/sanger_10"))),
                Some("sanger_5") => Ok(lines_from_bytes(include_bytes!("train/sanger_5"))),
                Some("start") => Ok(lines_from_bytes(include_bytes!("train/start"))),
                Some("start1") => Ok(lines_from_bytes(include_bytes!("train/start1"))),
                Some("stop") => Ok(lines_from_bytes(include_bytes!("train/stop"))),
                Some("stop1") => Ok(lines_from_bytes(include_bytes!("train/stop1"))),
                _ => Err(TrainingDataError::Io(filename.to_owned(), e)),
            }
        }
        Err(e) => Err(TrainingDataError::Io(filename.to_owned(), e)),
    }
}

fn next_line<R: BufRead>(
    filename: &PathBuf,
    lines: &mut Lines<R>,
) -> Result<String, TrainingDataError> {
    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)?
        .map_err(|e| TrainingDataError::Io(filename.to_path_buf(), e))
}

fn parse_float(
    filename: &PathBuf,
    part: &String,
    line: &String,
    string: &str,
) -> Result<f64, TrainingDataError> {
    f64::from_str(string).map_err(|e| {
        TrainingDataError::MalformedValue(filename.to_owned(), part.to_owned(), line.to_owned(), e)
    })
}

fn parse_float_col(
    filename: &PathBuf,
    part: &String,
    line: String,
    column: usize,
) -> Result<f64, TrainingDataError> {
    let v: Vec<&str> = line.split_whitespace().collect();
    parse_float(filename, part, &line, v[column])
}

fn read_transitions(global: &mut Global, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    let mut header;

    header = next_line(&filename, &mut lines)?;
    for _ in 0..NUM_TRANSITIONS {
        let line = next_line(&filename, &mut lines)?;
        let v = line.split_whitespace().collect::<Vec<&str>>();
        let value = parse_float(&filename, &header, &line, v[1])?.ln();
        match v[0] {
            "MM" => global.tr.mm = value,
            "MI" => global.tr.mi = value,
            "MD" => global.tr.md = value,
            "II" => global.tr.ii = value,
            "IM" => global.tr.im = value,
            "DD" => global.tr.dd = value,
            "DM" => global.tr.dm = value,
            "GE" => global.tr.ge = value,
            "GG" => global.tr.gg = value,
            "ER" => global.tr.er = value,
            "RS" => global.tr.rs = value,
            "RR" => global.tr.rr = value,
            "ES" => global.tr.es = value,
            "ES1" => global.tr.es1 = value,
            _ => Err(TrainingDataError::UnknownTransitionState)?,
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..ACGT {
        for j in 0..ACGT {
            global.tr_mi[i][j] =
                parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 2)?.ln();
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..ACGT {
        for j in 0..ACGT {
            global.tr_ii[i][j] =
                parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 2)?.ln();
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..State::COUNT {
        global.pi[i] =
            parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 1)?.ln();
    }

    Ok(())
}

fn read_m_transitions(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for p in 0..PERIOD {
            for c in 0..BI_ACGT {
                let line = next_line(&filename, &mut lines)?;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    locals[cg].e_m[p][c][e] = parse_float(&filename, &header, &line, v[e])?.ln();
                }
            }
        }
    }

    Ok(())
}

fn read_m1_transitions(
    locals: &mut Vec<Local>,
    filename: PathBuf,
) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for p in 0..PERIOD {
            for c in 0..BI_ACGT {
                let line = next_line(&filename, &mut lines)?;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    locals[cg].e_m1[p][c][e] = parse_float(&filename, &header, &line, v[e])?.ln();
                }
            }
        }
    }

    Ok(())
}

fn read_noncoding(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for e1 in 0..ACGT {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for e2 in 0..ACGT {
                locals[cg].tr_rr[e1][e2] = parse_float(&filename, &header, &line, v[e2])?.ln();
            }
        }
    }

    Ok(())
}

fn read_start(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_s[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_stop1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_e1[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_stop(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_e[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_start1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_s1[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_pwm(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    let mut line: String;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_s[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_e[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_s1[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_e1[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
    }

    Ok(())
}
