use sherloxome::cli::process_cli;

fn main() {
    env_logger::init();
    if let Err(e) = process_cli() {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
