use std::path::Path;
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-coverage"))
}

fn fixture() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/small.bam"))
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// (rname, numreads, covbases, coverage_pct, meandepth) per reference.
fn parse(s: &str) -> Vec<(String, u64, u64, f64, f64)> {
    s.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_string(),
                c[3].parse().unwrap(),
                c[4].parse().unwrap(),
                c[5].parse().unwrap(),
                c[6].parse().unwrap(),
            )
        })
        .collect()
}

// covbases is a union of CIGAR-aware aligned spans (correct under read
// overlap, soft-clips and indels); numbers must equal `samtools coverage`.
#[test]
fn coverage_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    let ours = bin().arg(fixture()).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    let theirs = Command::new("samtools")
        .arg("coverage")
        .arg(fixture())
        .output()
        .unwrap();
    assert!(theirs.status.success());

    let o = parse(&String::from_utf8_lossy(&ours.stdout));
    let t = parse(&String::from_utf8_lossy(&theirs.stdout));
    assert_eq!(o.len(), t.len(), "row count");
    for (a, b) in o.iter().zip(&t) {
        assert_eq!(a.0, b.0, "rname");
        assert_eq!(a.1, b.1, "numreads for {}", a.0);
        assert_eq!(a.2, b.2, "covbases for {}", a.0);
        assert!((a.3 - b.3).abs() < 1e-4, "coverage% {} vs {}", a.3, b.3);
        assert!((a.4 - b.4).abs() < 1e-4, "meandepth {} vs {}", a.4, b.4);
    }
}
