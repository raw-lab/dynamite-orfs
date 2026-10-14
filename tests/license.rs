//! DYNAMITE is CC BY-NC 4.0 throughout. An older Cargo.toml / README carried
//! GPL-3.0-or-later and kept coming back; these checks (and the compile-time
//! guard in `src/lib.rs`) fail if it does.

fn read(name: &str) -> String {
    std::fs::read_to_string(format!("{}/{}", env!("CARGO_MANIFEST_DIR"), name))
        .unwrap_or_else(|e| panic!("cannot read {name}: {e}"))
}

#[test]
fn cargo_toml_declares_cc_by_nc_4_0() {
    assert_eq!(env!("CARGO_PKG_LICENSE"), "CC-BY-NC-4.0");
    assert_eq!(dynamite::LICENSE, "CC-BY-NC-4.0");
    let toml = read("Cargo.toml");
    let lines: Vec<&str> = toml.lines().map(str::trim).filter(|l| l.starts_with("license")).collect();
    assert_eq!(lines, ["license = \"CC-BY-NC-4.0\""], "Cargo.toml licence lines: {lines:?}");
    assert!(!toml.contains("GPL"), "Cargo.toml mentions GPL");
}

#[test]
fn license_file_is_the_cc_by_nc_4_0_legal_code() {
    let text = read("LICENSE.txt");
    assert!(text.contains("Copyright (c) 2026 Richard Allen White III"));
    assert!(text.contains("SPDX identifier: CC-BY-NC-4.0"));
    assert!(text.contains("Creative Commons Attribution-NonCommercial 4.0 International Public License"));
    for section in 1..=8 {
        assert!(text.contains(&format!("Section {section} –")), "legal code section {section} missing");
    }
}

#[test]
fn readme_states_cc_by_nc_4_0() {
    let readme = read("README.md");
    assert!(readme.contains("[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](LICENSE.txt)"));
    assert!(!readme.contains("License-GPL"), "README carries a GPL licence badge");
    assert!(readme.contains("## 📄 License"));
    let footer = readme.lines().rev().find(|l| l.contains("<sub><b>DYNAMITE</b>")).expect("footer");
    assert!(footer.contains("CC BY-NC 4.0") && !footer.contains("GPL"), "footer: {footer}");
    let version_badge = format!("badge/version-{}-success", env!("CARGO_PKG_VERSION"));
    assert!(readme.contains(&version_badge), "README version badge is not {}", env!("CARGO_PKG_VERSION"));
}
