void HCalMapReorder(const char* inputfile, const char* outputfile, const double gain_multiplier = 1.0){
  ifstream infile(inputfile);
  ofstream outfile(outputfile);

  int idx_in, row, col, idx_out;
  double gain;
  double gain_array[288];
  idx_in = 0;
  while(infile.good() || idx_in<288){
    infile >> gain;
    col = idx_in%12;
    row = (idx_in-col)/12;
    cout << " idx_in " << idx_in << " row " << row << " col " << col;
    col = 12-1-col;
    idx_out = row*12 + col;
    cout << " => col " << col << " idx_out " << idx_out << " gain " << gain << endl;
    gain_array[idx_out] = gain*gain_multiplier;
    idx_in++;
  }
  
  for(idx_out = 0; idx_out<288; idx_out++){
    outfile << gain_array[idx_out] << " ";
  }
  outfile << endl;
}
