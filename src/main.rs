use std::fs::{File, write};
use std::io::{self, BufRead};
use std::path::Path;

fn main() {
    let data: Vec<String> = make_clean_fasta_data("rosalind_edta.txt"); //read in fasta data
    //println!("{:?}",data); //debug

    // get two strings to align
    let x: &String = &data[0];
    let y: &String = &data[1];

    // Needleman-Wunsch algorithm:
    let (z,w,distance) = alignment(x, y);
    let ans: String = format!("{}\n{}\n{}",distance,z,w);
    println!("{}",ans);
    write("ans.txt",ans).expect("Should write to file");
}

fn alignment(x: &String, y: &String) -> (String, String, usize) {
    // construct matrix
    let mut mat: Vec<Vec<usize>> = vec![vec![0usize;x.len()+1];y.len()+1];

    //preallocate first row
    for i in 0..x.len()+1 {
        mat[0][i] = i;
    }
    //preallocate first column
    for j in 0..y.len()+1 {
        mat[j][0] = j;
    }

    //iterate through matrix
    //current cell is i+1, j+1
    for j in 0..y.len() {
        for i in 0..x.len(){
            let cost:usize;
            // if it's a match
            if y[j..j+1] == x[i..i+1] {
                cost = 0;
            } else {
                cost = 1;
            }

            // score based on adjacent cells
            let top: usize = mat[j][i+1] + 1;
            let left: usize = mat[j+1][i] + 1;
            let diag: usize = mat[j][i] + cost;
            let adj_cells: [usize; 3] = [top,left,diag];

            // find max score
            let min_val: Option<&usize> = adj_cells.iter().min();
            match min_val {
                Some(cell) => mat[j+1][i+1] = *cell, //assign cell score
                None => ()
            }
        }
    }

    let edit_distance: usize = mat[y.len()][x.len()];

    /*
    for j in 0..mat.len(){ //debug, print matrix pretty
        println!("{:?}",mat[j])
    }
    */
    //find aligment strings
    //start at bottom-right of matrix, find path to top left
    //need to store candidates
    let mut z: String = String::new();
    let mut w: String = String::new();

    // CHOOSE PATHS
    // NEED WHILE TO DO (0,0) CELL
    let mut i: usize = x.len()-1;
    let mut j: usize = y.len()-1;
    while !(i == 0 && j == 0) {
        //current cell is [j+1][i+1]
        let top: usize = mat[j][i+1];
        let left: usize = mat[j+1][i];
        let diag: usize = mat[j][i];

        if x[i..i+1] == y[j..j+1] { //letters align
            z.push_str(&x[i..i+1]);
            w.push_str(&y[j..j+1]);
            i -= 1;
            j -= 1;
        } else if top < left && top <= diag { //top best path -> gap align to row letter
            z.push('-');
            w.push_str(&y[j..j+1]);
            j -= 1;
        } else if left < top && left <= diag { //left best path -> align gap to column letter
            w.push('-');
            z.push_str(&x[i..i+1]);
            i -= 1;
        } else { //diag best path ->  swap
            z.push_str(&x[i..i+1]);
            w.push_str(&y[j..j+1]);
            i -= 1;
            j -= 1;
        }
    }
    //push (0,0) cell
    z.push_str(&x[0..1]);
    w.push_str(&y[0..1]);
    z = z.chars().rev().collect();
    w = w.chars().rev().collect();


    return(z,w,edit_distance);
}

fn make_clean_fasta_data(filepath:&str) -> Vec<String> {
    let mut data:Vec<String> = vec![String::new();0]; // final output string vector
    let mut data_line:String = String::new(); // current strand

    // File hosts.txt must exist in the current path
    if let Ok(lines) = read_lines(filepath) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines.map_while(Result::ok) {
            //println!("{}",line); //debug
            let tag = &line[..1]; //FASTA ID line
            if tag == ">" { // if ID line
                data.push(data_line); // save previous strand //MAKES EMPTY FIRST SLOT
                data_line = String::new(); // start new strand
            } else {
                data_line.push_str(&line); //push file line to current strand
            }
        }
        data.push(data_line); // push final string to strand
    } else {
        println!("Can't read file!");
    }
    data.remove(0); // remove first empty slot
    //println!("{:?}",data); //debug
    return data;
}

// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}