use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_coverage(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-coverage");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    c.bench_function("rsomics-bam-coverage golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(bam.to_str().unwrap())
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bam_coverage);
criterion_main!(benches);
