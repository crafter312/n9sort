void readPIDzline(int quad, int crystal)
{
  int Nlines;
  ostringstream outstring;
  outstring << "/user/e25001/analysis/n9sort/zline/pid_" << quad << "_" << crystal  << ".zline";
  string sname = outstring.str();
  cout << sname << endl;
  ifstream zlines(sname.c_str());
  zlines >> Nlines;
  cout << Nlines << endl;
  TCutG *line[22];
  char name[10];
  int Z,A;

  for (int i=0;i<Nlines;i++)
  {
    int Npoints;
    zlines >> Z >> A;
    cout << "Z=" << Z << "  A=" << A << endl;
    zlines >> Npoints;
    cout << "finger #" << i << "  N=" << Npoints << endl;
    sprintf(name,"finger[%d]",i);
    line[i] = new TCutG(name,Npoints);

    for (int j=0;j<Npoints;j++)
    {
      double x,y;
      zlines >> x >> y;
      cout << "   " << x << " " << y << endl;
      line[i]->SetPoint(j,x,y);
    }
    line[i]->Draw();
  }
}
