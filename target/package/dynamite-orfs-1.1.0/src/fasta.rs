//! Minimal, streaming FASTA reader with transparent gzip support.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// A single FASTA record.
pub struct Record {
    /// The header id (first whitespace-delimited token after `>`).
    pub id: String,
    /// The full header line (without the leading `>`).
    pub desc: String,
    /// Raw sequence bytes (no newlines), as they appeared in the file.
    pub seq: Vec<u8>,
}

/// Read every record from a FASTA file. Files ending in `.gz` (or beginning
/// with the gzip magic bytes) are decompressed transparently.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Record>> {
    let path = path.as_ref();
    let file = File::open(path)
        .with_context(|| format!("could not open input file {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Peek the first two bytes to detect gzip regardless of extension. A
    // single `read` may legally return fewer bytes than asked for (pipes,
    // FIFOs), so keep reading until two bytes or EOF.
    let mut magic = [0u8; 2];
    let mut n = 0usize;
    while n < 2 {
        let got = reader.read(&mut magic[n..])?;
        if got == 0 {
            break;
        }
        n += got;
    }
    let is_gz = n == 2 && magic[0] == 0x1f && magic[1] == 0x8b;

    // Reconstruct a reader over the full stream (we consumed the magic bytes).
    let head: Vec<u8> = magic[..n].to_vec();
    let chained = io::Cursor::new(head).chain(reader);
    let boxed: Box<dyn Read> = if is_gz {
        Box::new(flate2::read::MultiGzDecoder::new(chained))
    } else {
        Box::new(chained)
    };

    parse(BufReader::new(boxed))
}

fn parse<R: BufRead>(reader: R) -> Result<Vec<Record>> {
    let mut records = Vec::new();
    let mut id = String::new();
    let mut desc = String::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut have = false;

    // Read raw byte lines: a header with a non-UTF-8 byte (e.g. Latin-1
    // strain names) used to abort the whole read with "invalid UTF-8".
    for line in reader.split(b'\n') {
        let bytes = line?;
        if bytes.first() == Some(&b'>') {
            if have {
                records.push(Record {
                    id: std::mem::take(&mut id),
                    desc: std::mem::take(&mut desc),
                    seq: std::mem::take(&mut seq),
                });
            }
            let text = String::from_utf8_lossy(&bytes[1..]);
            let header = text.trim_end();
            desc = header.to_string();
            id = header
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            have = true;
        } else if have {
            for &b in &bytes {
                if !b.is_ascii_whitespace() {
                    seq.push(b);
                }
            }
        }
        // lines before the first '>' are silently ignored
    }
    if have {
        records.push(Record { id, desc, seq });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_two_records() {
        let data = ">seq1 first\nACGT\nACGT\n>seq2\nTTTT\n";
        let recs = parse(BufReader::new(Cursor::new(data))).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "seq1");
        assert_eq!(recs[0].desc, "seq1 first");
        assert_eq!(recs[0].seq, b"ACGTACGT");
        assert_eq!(recs[1].id, "seq2");
        assert_eq!(recs[1].seq, b"TTTT");
    }

    #[test]
    fn non_utf8_header_and_crlf() {
        let data: &[u8] = b">s1 souche \xe9t\xe9\r\nACGT\r\nAC\n";
        let recs = parse(BufReader::new(Cursor::new(data))).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "s1");
        assert_eq!(recs[0].seq, b"ACGTAC");
    }
}
