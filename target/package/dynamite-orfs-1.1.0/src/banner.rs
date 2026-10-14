//! The DYNAMITE ASCII banner, shared by the CLI help and the `doctor` report.

/// Crate version, surfaced in the banner and `doctor`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The full banner: figlet wordmark with a lit stick-of-dynamite + mushroom-
/// cloud motif, plus a version / tagline footer.
pub fn banner() -> String {
    format!(
        "{}\nv{}  \u{00b7}  ORF / gene caller for all domains of life\nbacteria \u{00b7} archaea \u{00b7} phage \u{00b7} giant virus \u{00b7} crassphage \u{00b7} metagenome \u{00b7} eukaryote \u{00b7} reads\n",
        ART, VERSION
    )
}

/// Just the wordmark art (no footer): DYNAMITE figlet with a lit dynamite stick
/// erupting into a mushroom cloud.
pub const ART: &str = r#"                                                      _.-^^~~~~~^^-._
                                                  _.-~ ░░▒▒▓▓▓▓▒▒░░ ~-._
                                                .'  ░▒▓████████████▓▒░  '.
                                               ( ░▒▓██████████████████▓▒░ )
                                                '._ ░▒▓████████████▓▒░ _.'
                                                   '~-.░▒▓██████▓▒░.-~'
 ______   ___   _    _    __  __ ___ _____ _____      \ ░▒▓███▓▒░ /
|  _ \ \ / / \ | |  / \  |  \/  |_ _|_   _| ____|      \ ▒▓███▓▒ /
| | | \ V /|  \| | / _ \ | |\/| || |  | | |  _|   ╭──╮  \░▒▓█▓▒░/
| |_| || | | |\  |/ ___ \| |  | || |  | | | |___  ││▒▓██╼╾╼ ✸
|____/ |_| |_| \_/_/   \_\_|  |_|___| |_| |_____| ╰──╯"#;
