use edf::Reader;
use std::error::Error;

pub fn print_edf_signals(file_path: &str) -> Result<(), Box<dyn Error>> {
    let edf_file = Reader::from_path(file_path)?;
    println!("{}", edf_file.header);
    println!("Start datetime: {}", edf_file.header.start_datetime);
    println!("{}", edf_file.signal_header);
    println!("##data\n{}", edf_file.data);
    Ok(())
}
