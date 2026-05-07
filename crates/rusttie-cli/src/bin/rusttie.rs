use clap::Parser;

fn main() -> anyhow::Result<()> {
    rusttie_cli::run(rusttie_cli::Cli::parse())
}
