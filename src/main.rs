use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::path::Path;

mod csv_utils;
mod edf_utils;
mod models;
use csv_utils::read_ecg_data;
use models::EcgPoint;
use std::io::{self, Write};
// basic structure to hold our ECG data points

fn main() -> Result<(), Box<dyn Error>> {
    // getting the current directory
    let current_dir = std::env::current_dir()?;
    println!("Current directory: {:?}", current_dir);

    // file paths
    let input_path = current_dir.join("ecg.csv");
    let output_path = current_dir.join("positions.txt");

    println!("Reading from: {:?}", input_path);

    // reading the data
    let ecg_data = read_ecg_data(&input_path)?;

    if ecg_data.is_empty() {
        println!("No data found in the ECG file");
        return Ok(());
    }

    // detecting QRS complexes
    let qrs_positions = detect_qrs_complexes(&ecg_data);

    println!("Writing to: {:?}", output_path);
    println!("Found {} QRS complexes", qrs_positions.len());

    // writing results to file
    write_positions_to_file(&qrs_positions, &output_path)?;

    println!("Detection complete.");

    edf_utils::print_edf_signals("example.edf")?;

    Ok(())
}

fn detect_qrs_complexes(ecg_data: &[EcgPoint]) -> Vec<f64> {
    // if no data, return empty vector
    if ecg_data.is_empty() {
        return Vec::new();
    }

    // calculating sampling frequency
    let sample_period = if ecg_data.len() > 1 {
        ecg_data[1].time - ecg_data[0].time
    } else {
        0.005 // assuming 200Hz as default
    };
    let fs = 1.0 / sample_period;

    println!("Detected sampling frequency: {:.2} Hz", fs);

    // Process the data in segments to handle long ECGs
    let segment_size = (30.0 * fs) as usize; // 30 second segments
    let mut all_qrs_positions = Vec::new();

    let total_segments = (ecg_data.len() + segment_size - 1) / segment_size;

    for segment_idx in 0..total_segments {
        let start_idx = segment_idx * segment_size;
        let end_idx = std::cmp::min((segment_idx + 1) * segment_size, ecg_data.len());

        if end_idx - start_idx < (2.0 * fs) as usize {
            // Skip segments shorter than 2 seconds
            continue;
        }

        // Process this segment
        let segment_positions = process_segment(&ecg_data[start_idx..end_idx], fs);

        // Add segment positions to overall list
        for pos in segment_positions {
            all_qrs_positions.push(pos);
        }
    }

    // Sort all positions in case segments were processed out of order
    all_qrs_positions.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    // Remove duplicates from segment overlaps
    let mut final_positions = Vec::new();
    if !all_qrs_positions.is_empty() {
        final_positions.push(all_qrs_positions[0]);

        for i in 1..all_qrs_positions.len() {
            if all_qrs_positions[i] - final_positions.last().unwrap() >= 0.2 {
                // 200ms minimum
                final_positions.push(all_qrs_positions[i]);
            }
        }
    }

    final_positions
}

fn process_segment(segment_data: &[EcgPoint], fs: f64) -> Vec<f64> {
    // Extract voltage values
    let voltage: Vec<f64> = segment_data.iter().map(|point| point.voltage).collect();

    // Step 1: Normalization
    let mean: f64 = voltage.iter().sum::<f64>() / voltage.len() as f64;
    let normalized: Vec<f64> = voltage.iter().map(|&v| v - mean).collect();

    // Step 2: Find QRS complexes directly
    find_qrs_peaks_direct(&normalized, segment_data, fs)
}

fn find_qrs_peaks_direct(voltage: &[f64], ecg_data: &[EcgPoint], fs: f64) -> Vec<f64> {
    let mut qrs_positions = Vec::new();

    // Constants adjusted for physiological values
    let min_peak_distance = (0.5 * fs) as usize; // 500ms minimum between QRS complexes
    let window_size = (0.15 * fs) as usize; // 150ms search window

    // Calculate voltage variability
    let std_dev = calculate_std_dev(voltage);
    let threshold = 2.0 * std_dev; // Threshold based on signal variability

    // Find all potential peaks (both positive and negative)
    let mut peak_candidates = Vec::new();

    for i in window_size..(voltage.len() - window_size) {
        // Check if this point is a significant local extrema
        let is_positive_peak = voltage[i] > 0.0
            && voltage[i]
                >= *voltage[i - window_size..i]
                    .iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(&f64::NEG_INFINITY)
            && voltage[i]
                >= *voltage[i + 1..i + window_size]
                    .iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(&f64::NEG_INFINITY)
            && voltage[i].abs() > threshold;

        let is_negative_peak = voltage[i] < 0.0
            && voltage[i]
                <= *voltage[i - window_size..i]
                    .iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(&f64::INFINITY)
            && voltage[i]
                <= *voltage[i + 1..i + window_size]
                    .iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(&f64::INFINITY)
            && voltage[i].abs() > threshold;

        if is_positive_peak || is_negative_peak {
            peak_candidates.push((i, voltage[i].abs()));
        }
    }

    // Sort peaks by amplitude (largest first)
    peak_candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));

    // Filter peaks keeping only the strongest ones that are sufficiently far apart
    let mut selected_peaks = Vec::new();

    for &(idx, _) in &peak_candidates {
        // Check if this peak is far enough from all previously selected peaks
        let is_isolated = selected_peaks
            .iter()
            .all(|&prev_idx| (idx as isize - prev_idx as isize).abs() > min_peak_distance as isize);

        if is_isolated {
            selected_peaks.push(idx);
        }
    }

    // Sort selected peaks by position
    selected_peaks.sort_unstable();

    // For each selected peak, find the exact R or S wave
    // by looking in a window centered on the peak
    for &idx in &selected_peaks {
        // Determine the window to search for the R or S wave
        let search_window = (0.08 * fs) as usize; // 80ms window
        let start = if idx > search_window {
            idx - search_window
        } else {
            0
        };
        let end = std::cmp::min(idx + search_window, ecg_data.len());

        if start >= end {
            continue;
        }

        // Find the absolute maximum voltage in the window
        // This will be either the R peak or S peak
        let mut max_abs_idx = start;
        let mut max_abs_value = ecg_data[start].voltage.abs();

        for i in start + 1..end {
            let abs_voltage = ecg_data[i].voltage.abs();
            if abs_voltage > max_abs_value {
                max_abs_idx = i;
                max_abs_value = abs_voltage;
            }
        }

        qrs_positions.push(ecg_data[max_abs_idx].time);
    }

    qrs_positions
}

fn calculate_std_dev(data: &[f64]) -> f64 {
    let mean: f64 = data.iter().sum::<f64>() / data.len() as f64;
    let variance: f64 =
        data.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / data.len() as f64;
    variance.sqrt()
}

fn write_positions_to_file<P: AsRef<Path>>(positions: &[f64], path: P) -> io::Result<()> {
    let mut file = File::create(path)?;

    for &pos in positions {
        writeln!(file, "{:.6}", pos)?;
    }

    Ok(())
}
