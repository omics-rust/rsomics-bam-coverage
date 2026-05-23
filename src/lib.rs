#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use noodles::sam::alignment::record::cigar::op::Kind;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

/// Reads excluded from coverage by default, matching `samtools coverage`:
/// UNMAP (0x4) | SECONDARY (0x100) | QCFAIL (0x200) | DUP (0x400).
/// Supplementary (0x800) alignments are counted.
const DEFAULT_EXCL: u16 = 0x704;

#[derive(Debug, Default, Clone, Serialize)]
pub struct RefCoverage {
    pub name: String,
    pub length: u64,
    pub mapped_reads: u64,
    pub covered_bases: u64,
    pub mean_depth: f64,
    pub coverage_pct: f64,
}

pub fn compute_coverage(input: &Path, workers: NonZero<usize>) -> Result<Vec<RefCoverage>> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_seqs = header.reference_sequences();
    let ref_names: Vec<String> = ref_seqs.keys().map(ToString::to_string).collect();
    let ref_lens: Vec<u64> = ref_seqs
        .values()
        .map(|rs| rs.length().get() as u64)
        .collect();

    let mut per_ref_reads: HashMap<usize, u64> = HashMap::new();
    let mut per_ref_aligned: HashMap<usize, u64> = HashMap::new();
    let mut events: HashMap<usize, Vec<(usize, i64)>> = HashMap::new();

    for result in reader.records() {
        let record = result.map_err(RsomicsError::Io)?;
        if (record.flags().bits() & DEFAULT_EXCL) != 0 {
            continue;
        }

        let Some(tid) = record.reference_sequence_id().transpose().ok().flatten() else {
            continue;
        };
        let Some(start) = record
            .alignment_start()
            .transpose()
            .ok()
            .flatten()
            .map(|p| p.get())
        else {
            continue;
        };

        *per_ref_reads.entry(tid).or_insert(0) += 1;
        let chrom_events = events.entry(tid).or_default();
        let aligned = per_ref_aligned.entry(tid).or_insert(0);

        let mut ref_pos = start;
        for op in record.cigar().iter() {
            let op = op.map_err(RsomicsError::Io)?;
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    chrom_events.push((ref_pos, 1));
                    chrom_events.push((ref_pos + len, -1));
                    *aligned += len as u64;
                    ref_pos += len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len;
                }
                _ => {}
            }
        }
    }

    let mut result = Vec::with_capacity(ref_names.len());
    for (i, name) in ref_names.iter().enumerate() {
        let length = ref_lens[i];
        let mapped_reads = per_ref_reads.get(&i).copied().unwrap_or(0);
        let aligned = per_ref_aligned.get(&i).copied().unwrap_or(0);
        let covered_bases = events.get_mut(&i).map_or(0, |evs| union_length(evs));

        let mean_depth = if length > 0 {
            aligned as f64 / length as f64
        } else {
            0.0
        };
        let coverage_pct = if length > 0 {
            covered_bases as f64 / length as f64 * 100.0
        } else {
            0.0
        };
        result.push(RefCoverage {
            name: name.clone(),
            length,
            mapped_reads,
            covered_bases,
            mean_depth,
            coverage_pct,
        });
    }

    Ok(result)
}

/// Total number of distinct reference positions covered by at least one
/// interval (the union length), via a single sorted sweep of +1/-1 events.
/// After consuming all deltas at position `p`, `depth` is the coverage of the
/// half-open span `[p, next_event_pos)`; events balance to zero so `depth`
/// is 0 once the last event is consumed.
fn union_length(events: &mut [(usize, i64)]) -> u64 {
    events.sort_unstable();
    let mut covered: u64 = 0;
    let mut depth: i64 = 0;
    let mut i = 0;
    while i < events.len() {
        let p = events[i].0;
        while i < events.len() && events[i].0 == p {
            depth += events[i].1;
            i += 1;
        }
        if depth > 0 && i < events.len() {
            covered += (events[i].0 - p) as u64;
        }
    }
    covered
}

pub fn write_coverage(cov: &[RefCoverage], output: &mut dyn Write) -> Result<()> {
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    writeln!(
        out,
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth"
    )
    .map_err(RsomicsError::Io)?;
    for r in cov {
        writeln!(
            out,
            "{}\t1\t{}\t{}\t{}\t{:.4}\t{:.4}",
            r.name, r.length, r.mapped_reads, r.covered_bases, r.coverage_pct, r.mean_depth
        )
        .map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
