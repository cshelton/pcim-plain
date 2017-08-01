#include "traj.h"
#include <sstream>
#include <iomanip>

using namespace std;

void printstatic(ostream &os, const traj &tr) {
	for(int i=0;i<tr.sx.size();i++)
		os << "static " << i << ": " << tr.sx[i] << endl;
}

void printtime(ostream &os, double t, int colw, bool ishrs) {
	if (!ishrs) os << setw(colw) << t;
	else {
		int sec = floor(t*60*60);
		int min = sec/60;
		sec %= 60;
		int hr = min/60;
		min %= 60;
		int day = hr/24;
		hr %= 24;
		stringstream ss;
		ss << setfill('0');
		ss << day << "d " << setw(2) << hr << ':' << setw(2) << min;
		for(int i=ss.str().size();i<colw;i++) os << ' ';
		os << ss.str();
	}
}

void printdynamic(ostream &os, const traj &tr, bool incolumns, bool ishrs) {
	vector<decltype(tr[0].begin())> it;
	vector<int> varid;
	int ndone = 0;
	for(int i=0;i<tr.size();i++) {
		if (tr[i].empty()) continue;
		it.push_back(tr[i].begin());
		varid.push_back(i);
	}
	int colw = 10;
	os << "events:" << endl;
	auto osf = os.setf(ios_base::fixed);
	if (incolumns) {
		os << "time      ";
		for(int i=0;i<varid.size();i++)
			os << "|" << setw(colw) << varid[i];
		os << endl;
	}
	while(ndone<it.size()) {
		int i=-1;
		double t = numeric_limits<double>::infinity();
		for(int j=0;j<it.size();j++)
			if (it[j]!=tr[varid[j]].end() && it[j]->first<t) {
				t = it[j]->first;
				i = j;
			}
		if (incolumns) {
			printtime(os,t,colw,ishrs);
			for(int j=0;j<i;j++) os << "|          ";
			//os << "|" << setw(colw) << it[i]->second;
			os << "|" << setw(colw) << "X";
			for(int j=i+1;j<varid.size();j++) os << "|          ";
			os << endl;
		} else {
			os << setw(4) << i << " @ ";
			printtime(os,t,colw,ishrs);
			//os << " " << setw(colw) << it[i]->second;
			os << endl;
		}
		if (++it[i] == tr[varid[i]].end()) ndone++;
	}
	os.setf(osf);
}

void printtr(ostream &os, const traj &tr, bool incolumns, bool ishrs) {
	printstatic(os,tr);
	printdynamic(os,tr,incolumns,ishrs);
}

BOOST_CLASS_EXPORT_IMPLEMENT(vartraj)
BOOST_CLASS_EXPORT_IMPLEMENT(vartraj::rec)
BOOST_CLASS_EXPORT_IMPLEMENT(traj)
