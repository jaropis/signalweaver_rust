use crate::models::EcgPoint;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
pub fn read_ecg_data<P: AsRef<Path>>(path: P) -> Result<Vec<EcgPoint>, Box<dyn Error>> {
    // opening the file
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();
    let mut header_skipped = false;

    // reading each line
    for line in reader.lines() {
        let line = line?;

        // skipping the header
        if !header_skipped {
            header_skipped = true;
            continue;
        }

        // parsing each line
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() == 2 {
            let time = parts[0].trim().parse::<f64>()?;
            let voltage = parts[1].trim().parse::<f64>()?;

            data.push(EcgPoint { time, voltage });
        }
    }

    // Print total data points
    if !data.is_empty() {
        println!("Total data points: {}", data.len());
    }

    Ok(data)
}
