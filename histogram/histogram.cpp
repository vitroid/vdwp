/*
 * A simpler implementation of calcularting the f-value.  It just calculate the histogram of the potential energy.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <queue>
#include <cmath>
#include <memory>  // Add this include at the top

// Memory usage is too large.
// Doubt some memory leakage.

using namespace std;

typedef vector<double> Coord;
typedef vector<int> Point; // Point should be 3-digits; not 6-digits including orientation.
typedef map<Point, double> Mark;
typedef queue<Point> Queue;
//格子間隔。辞書のkeyは生の実数はまずいので整数3つで指定する。
const double intv = 0.20; // angstrom; 0.20 is enough for 3 digits
const double dv = intv * intv * intv;
const int diva = 40;              // must be 40
const double intva = M_PI / diva; //(180 / diva) degree

/*

 vmin         vmin+binwidth     vmin+2binwidth
   |                 |                 |                 |                 |                 |
   |.................|.................|.................|.................|.................|
         binwidth
          bin 0               1                 2              ....               nbins-1
*/

// automatic histogram recorder
class cHistogram
{
private:
  int nbins, filledmin, filledmax, lot;
  double vmin, binwidth, totalweight;
  vector<double> histo;

public:
  cHistogram(double binwidth_) : binwidth(binwidth_), totalweight(0.0)
  {
    histo = vector<double>(1, 0.0);
    lot = 10;
  }

  // Rule of Five implementation
  cHistogram(const cHistogram&) = default;
  cHistogram& operator=(const cHistogram&) = default;
  cHistogram(cHistogram&&) = default;
  cHistogram& operator=(cHistogram&&) = default;
  ~cHistogram() = default;

  void accum(double value, double weight)
  {
    if (totalweight == 0)
    {
      // first process
      int bin = static_cast<int>(std::floor(value / binwidth));
      nbins = 1;
      vmin = bin * binwidth;
      filledmin = filledmax = 0;
    }
    if (value < vmin)
    {
      // how many bins should be added?
      int bin = (int)floor((vmin - value) / binwidth) * 2;
      while (lot < bin)
      {
        lot *= 2;
      }
      bin = lot;
      histo.resize(nbins + bin);
      for (int i = nbins - 1; i >= 0; i--)
      {
        histo[i + bin] = histo[i];
      }
      for (int i = 0; i < bin; i++)
      {
        histo[i] = 0.0;
      }
      nbins += bin;
      filledmin += bin;
      filledmax += bin;
      cout << "RESIZE -" << bin << endl;
      vmin -= bin * binwidth;
    }
    if (vmin + nbins * binwidth < value)
    {
      int bin = (int)floor((value - (vmin + nbins * binwidth)) / binwidth) * 2 + 1;
      while (lot < bin)
      {
        lot *= 2;
      }
      bin = lot;
      histo.resize(nbins + bin);
      nbins += bin;
      cout << "RESIZE +" << bin << endl;
    }
    // int is a mathematically ugry function.
    int bin = (int)floor((value - vmin) / binwidth);
    histo[bin] += weight;
    totalweight += weight;
    if (bin < filledmin)
    {
      filledmin = bin;
    }
    if (filledmax < bin)
    {
      filledmax = bin;
    }
  }
  void dump(ostream &fout, const string& tag)
  {
    if (totalweight == 0)
      return;
    if (tag == "@HIST")
    {
      cout << tag << endl;
      cout << 1 << endl; // dimension is always 1
      cout << filledmax - filledmin + 1 << " " << vmin + binwidth * filledmin << " " << binwidth << endl;
      for (int i = filledmin; i <= filledmax; i++)
      {
        cout << histo[i] << endl;
      }
    }
    if (tag == "@SHST")
    {
      // sparse histogram format
      cout << tag << endl;
      cout << 1 << endl; // dimension is always 1
      cout << vmin + binwidth * filledmin << " " << binwidth << endl;
      for (int i = filledmin; i <= filledmax; i++)
      {
        if (histo[i] != 0.0)
        {
          cout << i - filledmin << " " << histo[i] << endl;
        }
      }
      cout << "-1 0.0" << endl; // terminator
    }
  }
};

typedef struct
{
  double dEps;
  double dSig;
  double dCharge;
} sInteraction;

class cMolecule
{
private:
public:
  vector<sInteraction> intr; // nSite
  int nSite;
  vector<double> dMass;  // nSite
  vector<double> dCoord; // intramolecular coord (in the present case, Z coord only)
  cMolecule(istream &fin, string tag)
  {
    string buf;
    // fgets( buf, sizeof(buf), file ); // ID08 label
    getline(fin, buf);
    nSite = atoi(buf.c_str());
    intr = vector<sInteraction>(nSite);
    dMass = vector<double>(nSite);
    dCoord = vector<double>(nSite * 3);
    for (int i = 0; i < nSite; i++)
    {
      getline(fin, buf);
      double x, y, z;
      if (tag == "@DEFR")
        sscanf(buf.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &dMass[i]);
      else
      {
        x = y = z = 0.0;
        sscanf(buf.c_str(), "%lf", &dMass[i]);
      }
      dCoord[i * 3 + 0] = x;
      dCoord[i * 3 + 1] = y;
      dCoord[i * 3 + 2] = z;
    }
    for (int i = 0; i < nSite; i++)
    {
      getline(fin, buf);
      double eps, sig, ch;
      sscanf(buf.c_str(), "%lf %lf %lf", &eps, &sig, &ch);
      intr[i].dEps = eps;
      intr[i].dSig = sig;
      intr[i].dCharge = ch;
    }
  }
};

// Physical propertires
const double cdNA = 6.022045e23L;    // 1
const double cde = 1.602176565e-19L; // Coulomb
const double cdeps0 = 8.854e-12;     // permittivity
const double cdh = 6.626176e-34L;    // J s

vector<double> load_NX4A(istream &fin,
                         const cMolecule &lattice,
                         const vector<double> box)
{
  string buf;
  getline(fin, buf);
  int nmol = atoi(buf.c_str());
  vector<double> coord(nmol * lattice.nSite * 3);
  for (int mol = 0; mol < nmol; mol++)
  {
    getline(fin, buf);
    Coord pos(3);
    double a, b, c, d;
    sscanf(buf.c_str(), "%lf %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &a, &b, &c, &d);
    if (box.size() == 3)
    {
      for (int dim = 0; dim < 3; dim++)
      {
        pos[dim] -= rint(pos[dim] / box[dim]) * box[dim];
      }
    }
    double sp11 = (a * a + b * b - (c * c + d * d));
    double sp12 = -2.0 * (a * d + b * c);
    double sp13 = 2.0 * (b * d - a * c);
    double sp21 = 2.0 * (a * d - b * c);
    double sp22 = a * a + c * c - (b * b + d * d);
    double sp23 = -2.0 * (a * b + c * d);
    double sp31 = 2.0 * (a * c + b * d);
    double sp32 = 2.0 * (a * b - c * d);
    double sp33 = a * a + d * d - (b * b + c * c);
    for (int site = 0; site < lattice.nSite; site++)
    {
      double x = lattice.dCoord[site * 3 + 0];
      double y = lattice.dCoord[site * 3 + 1];
      double z = lattice.dCoord[site * 3 + 2];
      coord[(site + mol * lattice.nSite) * 3 + 0] = pos[0] + sp11 * x + sp12 * y + sp13 * z;
      coord[(site + mol * lattice.nSite) * 3 + 1] = pos[1] + sp21 * x + sp22 * y + sp23 * z;
      coord[(site + mol * lattice.nSite) * 3 + 2] = pos[2] + sp31 * x + sp32 * y + sp33 * z;
    }
  }
  return coord;
}

vector<double> load_AR3A(istream &fin,
                         const cMolecule &lattice,
                         const vector<double> box)
{
  string buf;
  getline(fin, buf);
  int nmol = atoi(buf.c_str());
  vector<double> coord(nmol * lattice.nSite * 3);
  for (int mol = 0; mol < nmol; mol++)
  {
    getline(fin, buf);
    Coord pos(3);
    sscanf(buf.c_str(), "%lf %lf %lf", &pos[0], &pos[1], &pos[2]);
    if (box.size() == 3)
    {
      for (int dim = 0; dim < 3; dim++)
      {
        pos[dim] -= rint(pos[dim] / box[dim]) * box[dim];
      }
    }
    coord[(0 + mol * lattice.nSite) * 3 + 0] = pos[0];
    coord[(0 + mol * lattice.nSite) * 3 + 1] = pos[1];
    coord[(0 + mol * lattice.nSite) * 3 + 2] = pos[2];
  }
  return coord;
}

// interaction from an atom of the guest to lattice molecules
// energy in kJ/mol
double
interaction(int nlattice, const vector<double> &latticeSites, const Coord &guest, const cMolecule &lattice, const sInteraction &intr)
{
  double ep = 0.0;
  int nlatticesite = lattice.nSite;
  for (int site = 0; site < nlatticesite; site++)
  {
    double eps = sqrt(lattice.intr[site].dEps * intr.dEps);
    double sig = (lattice.intr[site].dSig + intr.dSig) / 2;
    double chprod = lattice.intr[site].dCharge * intr.dCharge;
    if ((chprod != 0) || (eps != 0.0))
    {
      for (int mol = 0; mol < nlattice; mol++)
      {
        double d = 0.0;
        for (int dim = 0; dim < 3; dim++)
        {
          double r = latticeSites[(site + mol * lattice.nSite) * 3 + dim] - guest[dim];
          d += r * r;
        }
        double e1 = 0;
        if (eps != 0.0)
        {
          double di = (sig * sig) / d;
          double d6 = di * di * di;
          double d12 = d6 * d6;
          e1 = 4.0 * eps * (d12 - d6);
        }
        double e2 = 0;
        if (chprod != 0.0)
        {
          double r = sqrt(d);
          e2 = chprod * cde * cde / (4. * M_PI * cdeps0 * r * 1e-10) * cdNA * 1e-3; // kJ/mol
        }
        ep += e1 + e2;
      }
    }
  }
  return ep;
}

// energy in kJ/mol
// assume all the atoms in a molecule are on Z axis.
double
interaction2_rodlike(int nlattice, const vector<double> &latticeSites, const Point &point, const double theta, const double phi, const cMolecule &lattice, const cMolecule &molec)
{
  double ep;
  Coord axis(3, 0.0);
  // Computer Simulation of Liquids p.133
  axis[0] = cos(phi) * sin(theta);
  axis[1] = sin(phi) * sin(theta);
  axis[2] = cos(theta);
  Coord atom(3, 0.0);
  ep = 0.0;
  // int nlatticesite = lattice.nSite;
  for (int site = 0; site < molec.nSite; site++)
  {
    // assume atoms are on Z axis
    for (int i = 0; i < 3; i++)
    {
      atom[i] = point[i] * intv + axis[i] * molec.dCoord[site * 3 + 2];
    }
    ep += interaction(nlattice, latticeSites, atom, lattice, molec.intr[site]);
  }
  return ep;
}

// energy in kJ/mol
//田中先生の示唆により回転行列を変更。
double
interaction2_rigidbody(int nlattice, const vector<double> &latticeSites, const Point &point, double theta, double phi, double psi, const cMolecule &lattice, const cMolecule &molec)
{
  // Computer Simulation of Liquids p.133
  double sinth = sin(theta);
  double costh = cos(theta);
  double sinph = sin(phi);
  double cosph = cos(phi);
  double sinps = sin(psi);
  double cosps = cos(psi);
  // Goldstein p.
  double a11 = cosps * cosph - costh * sinph * sinps;
  double a12 = -cosps * sinph - costh * cosph * sinps;
  double a13 = sinps * sinth;
  double a21 = +sinps * cosph + costh * sinph * cosps;
  double a22 = -sinps * sinph + costh * cosph * cosps;
  double a23 = -cosps * sinth;
  double a31 = sinth * sinph;
  double a32 = +sinth * cosph;
  double a33 = costh;
  Coord atom(3, 0.0);
  double ep = 0.0;
  for (int site = 0; site < molec.nSite; site++)
  {
    // assume atoms are on Z axis
    atom[0] = point[0] * intv + a11 * molec.dCoord[site * 3 + 0] + a12 * molec.dCoord[site * 3 + 1] + a13 * molec.dCoord[site * 3 + 2];
    atom[1] = point[1] * intv + a21 * molec.dCoord[site * 3 + 0] + a22 * molec.dCoord[site * 3 + 1] + a23 * molec.dCoord[site * 3 + 2];
    atom[2] = point[2] * intv + a31 * molec.dCoord[site * 3 + 0] + a32 * molec.dCoord[site * 3 + 1] + a33 * molec.dCoord[site * 3 + 2];
    ep += interaction(nlattice, latticeSites, atom, lattice, molec.intr[site]);
  }
  return ep;
}

double
interaction2_rigidbody_original_euler(int nlattice, const vector<double> &latticeSites, const Point &point, const cMolecule &lattice, const cMolecule &molec)
{
  double theta = point[3] * intva;
  double phi = point[4] * intva;
  double psi = point[5] * intva;
  // Computer Simulation of Liquids p.133
  double sinth = sin(theta);
  double costh = cos(theta);
  double sinph = sin(phi);
  double cosph = cos(phi);
  double sinps = sin(psi);
  double cosps = cos(psi);
  double a11 = cosps * cosph - costh * sinph * sinps;
  double a12 = cosps * sinph + costh * cosph * sinps;
  double a13 = sinps * sinth;
  double a21 = -sinps * cosph - costh * sinph * cosps;
  double a22 = -sinps * sinph + costh * cosph * cosps;
  double a23 = cosps * sinth;
  double a31 = sinth * sinph;
  double a32 = -sinth * cosph;
  double a33 = costh;
  Coord atom(3, 0.0);
  double ep = 0.0;
  for (int site = 0; site < molec.nSite; site++)
  {
    // assume atoms are on Z axis
    atom[0] = point[0] * intv + a11 * molec.dCoord[site * 3 + 0] + a12 * molec.dCoord[site * 3 + 1] + a13 * molec.dCoord[site * 3 + 2];
    atom[1] = point[1] * intv + a21 * molec.dCoord[site * 3 + 0] + a22 * molec.dCoord[site * 3 + 1] + a23 * molec.dCoord[site * 3 + 2];
    atom[2] = point[2] * intv + a31 * molec.dCoord[site * 3 + 0] + a32 * molec.dCoord[site * 3 + 1] + a33 * molec.dCoord[site * 3 + 2];
    ep += interaction(nlattice, latticeSites, atom, lattice, molec.intr[site]);
  }
  return ep;
}

typedef vector<string> VS;
VS split(string s, string c)
{
  VS ret;
  for (int i = 0, n; i <= s.length(); i = n + 1)
  {

    n = s.find_first_of(c, i);
    if (n == string::npos)
      n = s.length();
    string tmp = s.substr(i, n - i);
    ret.push_back(tmp);
  }
  return ret;
}

void probe(Queue &q,
           int nlattice,
           const vector<double> &latticeSites,
           Mark &mark,
           const cMolecule &lattice,
           const cMolecule &molec,
           int dimen,
           cHistogram &histo)
{
  Point point = q.front();
  q.pop();
  // return if the grid point is already calculated.
  if (mark.find(point) != mark.end())
  {
    return;
  }
  double emin = 0.0;
  if (dimen == 0)
  {
    // monatomic guest
    vector<double> atom(3, 0.0);
    for (int dim = 0; dim < 3; dim++)
    {
      atom[dim] = point[dim] * intv;
    }
    double ep = 0;
    for (int site = 0; site < molec.nSite; site++)
    {
      ep += interaction(nlattice, latticeSites, atom, lattice, molec.intr[site]);
    }
    emin = ep;
    if (ep < 0.0)
      histo.accum(ep, dv / 1e30);
    mark[point] = 1;
  }
  else if (dimen == 1)
  {
    // rod-like
    for (int itheta = 0; itheta < diva; itheta++)
    { // 180 degree in theta
      double theta = itheta * intva;
      double weight = sin(theta); // OK
      for (int iphi = 0; iphi < diva * 2; iphi++)
      { // 360 degree in phi
        double phi = iphi * intva;
        double energy = interaction2_rodlike(nlattice, latticeSites, point, theta, phi, lattice, molec);
        if (energy < emin)
        {
          emin = energy;
        }
        if (energy < 0.0)
          histo.accum(energy, weight * intva * intva * dv / 1e30);
      }
    }
    mark[point] = 1;
  }
  else if (dimen == 3)
  {
    // rigid body
    for (int itheta = 0; itheta < diva; itheta++)
    { // 180 degree in theta
      double theta = itheta * intva;
      double weight = sin(theta); // OK
      for (int iphi = 0; iphi < diva * 2; iphi++)
      { // 360 degree in phi
        double phi = iphi * intva;
        for (int ipsi = 0; ipsi < diva * 2; ipsi++)
        { // 360 degree in psi
          double psi = ipsi * intva;
          double energy = interaction2_rigidbody(nlattice, latticeSites, point, theta, phi, psi, lattice, molec);
          if (energy < emin)
          {
            emin = energy;
          }
          if (energy < 0.0)
            histo.accum(energy, weight * intva * intva * intva * dv / 1e30);
        }
      }
    }
    mark[point] = 1;
  }
  // if the energy exceeds the threshold, terminate calculation.
  if (emin < 0)
  {
    // Otherwise, scan the lattice recursively.
    Point p;
    p = point;
    p[0] += 1;
    q.push(p);
    p = point;
    p[0] -= 1;
    q.push(p);
    p = point;
    p[1] += 1;
    q.push(p);
    p = point;
    p[1] -= 1;
    q.push(p);
    p = point;
    p[2] += 1;
    q.push(p);
    p = point;
    p[2] -= 1;
    q.push(p);
  }
}

int MolecularShape(const cMolecule &molec)
{
  vector<double> moment_of_inertia(3, 0.0);
  // double sym = 0.0;
  double mass = 0.0;
  for (int site = 0; site < molec.nSite; site++)
  {
    mass += molec.dMass[site];
    moment_of_inertia[0] += molec.dMass[site] * (molec.dCoord[site * 3 + 1] * molec.dCoord[site * 3 + 1] + molec.dCoord[site * 3 + 2] * molec.dCoord[site * 3 + 2]);
    moment_of_inertia[1] += molec.dMass[site] * (molec.dCoord[site * 3 + 2] * molec.dCoord[site * 3 + 2] + molec.dCoord[site * 3 + 0] * molec.dCoord[site * 3 + 0]);
    moment_of_inertia[2] += molec.dMass[site] * (molec.dCoord[site * 3 + 0] * molec.dCoord[site * 3 + 0] + molec.dCoord[site * 3 + 1] * molec.dCoord[site * 3 + 1]);
    // sym += molec.dMass[site] * (molec.dCoord[site * 3 + 2] * molec.dCoord[site * 3 + 2] * molec.dCoord[site * 3 + 2]);
  }
  // cout << sym << " SYM" << endl;
  int molecular_shape = 0; // molecular shape
  if (moment_of_inertia[2] == 0.0)
  {
    if (moment_of_inertia[1] > 0.0)
    {
      molecular_shape = 1;
    }
  }
  else
  {
    molecular_shape = 3;
  }
  return molecular_shape;
}

void histogram(cHistogram &histo,
               int nlattice,
               const vector<double> &latticeSites,
               const cMolecule &lattice,
               const cMolecule &molec)
{
  cout << molec.nSite << endl;

  auto dimen = MolecularShape(molec);

  cout << dimen << " DIMEN" << endl;
  // sym should be zero if symmetric
  Mark mark;
  Point origin(3, 0); // three dimensional location
  Queue q;
  q.push(origin);
  int count = 0;
  int intv = 1;
  while (!q.empty())
  {
    probe(q, nlattice, latticeSites, mark, lattice, molec, dimen, histo);
    count += 1;

    if (intv == count)
    {
      int qlen = q.size();
      cout << count << "(" << qlen << "):" << endl;
      intv *= 2;
    }
    if (count > 3000000)
    {
      cerr << "TOO MANY SAMPLES" << count << endl;
      return;
    }
  }
  typedef Mark::const_iterator CI;
  for (CI p = mark.begin(); p != mark.end(); ++p)
  {
    Point v = p->first;
    cout << v[0] << "," << v[1] << "," << v[2] << ":" << p->second << endl;
  }
}

int countCoord(int nlattice,
               const vector<double> &latticeSites,
               const cMolecule &lattice,
               double radius)
{
  int nnei = 0;
  for (int mol = 0; mol < nlattice; mol++)
  {
    double mass = 0.0;
    vector<double> pos(3, 0.0);
    for (int site = 0; site < lattice.nSite; site++)
    {
      mass += lattice.dMass[site];
      for (int dim = 0; dim < 3; dim++)
      {
        double r = latticeSites[(site + mol * lattice.nSite) * 3 + dim];
        pos[dim] += lattice.dMass[site] * r;
      }
    }
    double dd = 0.0;
    for (int dim = 0; dim < 3; dim++)
    {
      pos[dim] /= mass;
      double r = pos[dim];
      dd += r * r;
    }
    if (dd < radius * radius)
    {
      nnei++;
    }
  }
  return nnei;
}

// interaction parameters
const double cdJ2K = 0.120273; // K/J

// const double cdEpsTIP4Ptk = 0.6487; // kJ/mol
// const double cdSigTIP4Ptk = 3.154L; //Angstrom

#include <opt.h>

int main(int argc, char *argv[])
{
  cout.precision(17);
  char *guestID = "LJME____";
  char* latticeID = "TIP4P   ";
  map<string, std::unique_ptr<cMolecule>> defr;  // Changed to unique_ptr
  OptRegister(&guestID, 'i', "guest", "8-letter guest ID");
  OptRegister(&latticeID, 'w', "water", "8-letter water ID");
  opt(&argc, &argv);

  cHistogram histo(0.01);
  vector<double> box;
  cout << "# cage size; f-value/(kJ/mol)" << endl;
  string tag;
  while (getline(cin, tag))
  {
    if (tag.size() > 4 && ((tag.substr(0, 5) == "@DEFR") || (tag.substr(0, 5) == "@DEFP")))
    {
      string buf;
      getline(cin, buf);
      buf += "        ";
      buf = buf.substr(0, 8);
      cout << tag << ":" << buf << endl;
      defr[buf] = std::make_unique<cMolecule>(cin, tag.substr(0, 5));  // Using make_unique
      cout << buf << "//" << defr[buf].get() << endl;
    }
    // else if (tag.size() > 4 && (tag.substr(0, 5) == "@ID08"))
    // {
    //   string buf;
    //   getline(cin, buf);
    //   buf += "        ";
    //   latticeID = buf.substr(0, 8);
    // }
    else if (tag.size() > 4 && (tag.substr(0, 5) == "@BOX3"))
    {
      string buf;
      getline(cin, buf);
      box = vector<double>(3, 0.0);
      sscanf(buf.c_str(), "%lf %lf %lf", &box[0], &box[1], &box[2]);
    }
    else if (tag.size() > 4 && (tag.substr(0, 5) == "@NX4A" ||
                                tag.substr(0, 5) == "@AR3A"))
    {
      // read lattice coord
      cMolecule* lattice = defr[latticeID].get();  // Get raw pointer when needed
      cout << latticeID << endl;
      cout << lattice << endl;
      vector<double> latticeSites;
      if (tag.substr(0, 5) == "@NX4A")
        latticeSites = load_NX4A(cin, *lattice, box);
      else
        latticeSites = load_AR3A(cin, *lattice, box);
      int nlattice = latticeSites.size() / (3 * lattice->nSite);
      int nnei = countCoord(nlattice, latticeSites, *lattice, 6.0);
      cout << nnei << " cagesize" << endl;
      string molec(guestID);
      cout << guestID << endl;
      cout << molec << endl;
      cout << defr[molec].get() << endl;
      histogram(histo, nlattice, latticeSites, *lattice, *defr[molec]);
    }
  }
  histo.dump(cout, "@SHST");
}
// 1st test: TIP4PLJ in TIP4PLJ...OK
// LJ in TIP4PLJ .. OK
// propane (3sites) in TIP4PLJ ..OK
