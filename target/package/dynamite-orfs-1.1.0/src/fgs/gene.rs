extern crate thiserror;
use thiserror::Error;

use crate::fgs::dna::Nuc::{A, C, G, T};
use crate::fgs::dna::{trinucleotide, Nuc, ANTI_CODON_CODE, CODON_CODE};

#[derive(Debug, PartialEq)]
pub struct ReadPrediction {
    pub head: Vec<u8>,
    pub genes: Vec<Gene>,
}

impl ReadPrediction {
    pub fn new(head: Vec<u8>) -> Self {
        ReadPrediction {
            head: head,
            genes: vec![],
        }
    }

    pub fn append_to(
        &self,
        aabuf: &mut Option<Vec<u8>>,
        metabuf: &mut Option<Vec<u8>>,
        dnabuf: &mut Option<Vec<u8>>,
        formatted: bool,
        whole_genome: bool,
    ) -> Result<(), GeneError> {
        if let Some(metabuf) = metabuf {
            self.meta(&mut *metabuf)?;
        }
        if let Some(dnabuf) = dnabuf {
            self.dna(&mut *dnabuf, formatted)?;
        }
        if let Some(aabuf) = aabuf {
            self.protein(&mut *aabuf, whole_genome)?;
        }
        Ok(())
    }

    pub fn meta(&self, buf: &mut Vec<u8>) -> Result<(), GeneError> {
        if !self.genes.is_empty() {
            buf.append(&mut format!(">{}\n", std::str::from_utf8(&self.head)?).into_bytes())
        }
        for gene in &self.genes {
            gene.meta(buf);
        }
        Ok(())
    }

    pub fn gff(&self, buf: &mut Vec<u8>) -> Result<(), GeneError> {
        if !self.genes.is_empty() {
            let head = std::str::from_utf8(&self.head)?;
            for gene in &self.genes {
                gene.gff(buf, &head);
            }
        }
        Ok(())
    }

    pub fn dna(&self, buf: &mut Vec<u8>, formatted: bool) -> Result<(), GeneError> {
        for gene in &self.genes {
            gene.dna(buf, &self.head, formatted)?;
        }
        Ok(())
    }

    pub fn protein(&self, buf: &mut Vec<u8>, whole_genome: bool) -> Result<(), GeneError> {
        for gene in &self.genes {
            gene.protein(buf, &self.head, whole_genome)?;
        }
        Ok(())
    }
}

#[derive(Debug, PartialEq)]
pub struct Gene {
    pub start: usize,
    pub end: usize,
    pub frame: usize,
    pub score: f64,
    pub dna: Vec<Nuc>,
    pub forward_strand: bool,
    pub inserted: Vec<usize>,
    pub deleted: Vec<usize>,
}

impl Gene {
    pub fn meta(&self, buf: &mut Vec<u8>) {
        buf.append(
            &mut format!(
                "{}\t{}\t{}\t{}\t{:.6}\tI:{}\tD:{}\n",
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                self.frame,
                self.score,
                self.inserted
                    .iter()
                    .map(|i: &usize| { format!("{},", i) })
                    .collect::<String>(),
                self.deleted
                    .iter()
                    .map(|i: &usize| { format!("{},", i) })
                    .collect::<String>()
            )
            .into_bytes(),
        );
    }

    pub fn gff(&self, buf: &mut Vec<u8>, head: &str) {
        buf.append(
            &mut format!(
                "{}\tFGS\tCDS\t{}\t{}\t.\t{}\t{}\tID={}_{}_{}_{};product=predicted protein\n",
                head,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                self.frame - 1,
                head,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' }
            )
            .into_bytes(),
        );
    }

    pub fn dna(&self, buf: &mut Vec<u8>, head: &Vec<u8>, formatted: bool) -> Result<(), GeneError> {
        let dna: Vec<u8> = match (self.forward_strand, formatted) {
            (true, true) => self.dna.iter().map(|&n| u8::from(n)).collect(),
            (true, false) => self
                .dna
                .iter()
                .filter(|n| !n.is_insertion())
                .map(|&n| u8::from(n))
                .collect(),
            (false, true) => self.dna.iter().rev().map(|&n| u8::from(n.rc())).collect(),
            (false, false) => self
                .dna
                .iter()
                .rev()
                .filter(|n| !n.is_insertion())
                .map(|&n| u8::from(n.rc()))
                .collect(),
        };

        buf.append(
            &mut format!(
                ">{}_{}_{}_{}\n{}\n",
                std::str::from_utf8(head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                std::str::from_utf8(&dna)?,
            )
            .into_bytes(),
        );

        Ok(())
    }

    pub fn protein(
        &self,
        buf: &mut Vec<u8>,
        head: &Vec<u8>,
        whole_genome: bool,
    ) -> Result<(), GeneError> {
        let dna = self
            .dna
            .iter()
            .filter(|n| !n.is_insertion())
            .map(|&n| n)
            .collect::<Vec<Nuc>>();
        let mut protein: Vec<u8> = if self.forward_strand {
            dna.chunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        } else {
            dna.rchunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| ANTI_CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        };
        if protein.last() == Some(&b'*') {
            protein.pop();
        }

        // alternative start codons still encode for Met
        // E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG and one or two others (e.g., an AUU and possibly a CUG)
        // only consider two major alternative ones, GTG and TTG
        if whole_genome {
            if self.forward_strand {
                let s = trinucleotide(&self.dna);
                if s == trinucleotide(&[G, T, G]) || s == trinucleotide(&[T, T, G]) {
                    protein[0] = b'M';
                }
            } else {
                let s = trinucleotide(self.dna.get(self.dna.len() - 3..).unwrap());
                if s == trinucleotide(&[C, A, C]) || s == trinucleotide(&[C, A, A]) {
                    protein[0] = b'M';
                }
            }
        }

        buf.append(
            &mut format!(
                ">{}_{}_{}_{}\n{}\n",
                std::str::from_utf8(head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                std::str::from_utf8(&protein)?,
            )
            .into_bytes(),
        );
        Ok(())
    }
}

#[derive(Error, Debug)]
pub enum GeneError {
    #[error("could not convert header back to UTF-8")]
    Utf8Error(#[from] std::str::Utf8Error),
}
