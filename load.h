#ifndef LOAD_H
#define LOAD_H

#include <string>
#include <fstream>
#include <sstream>
#include "traj.h"
#include <utility>
#include <cmath>
#include "serial.h"
#include <boost/regex.hpp>

inline std::vector<std::string> strsplit(const std::string &s, char delim) {
	std::vector<std::string> ret;
	std::istringstream ss(s);
	for(std::string tok; std::getline(ss,tok,delim); )
		ret.emplace_back(tok);
	return ret;
}

inline void akiloadep(std::string fn, traj &tr, int nsvar, int ndvar,
			int wtvarnum, const std::set<int> &tonormbywt,
			const std::map<int,int> *smap=nullptr,
			const std::map<int,int> *dmap=nullptr) {
	std::ifstream s(fn.c_str());
	if (!s.good()) return; // error out?
	double maxtime = 0.0;
	tr.resize(ndvar);
	tr.sx.resize(nsvar,std::numeric_limits<double>::quiet_NaN());
	double wt = -1.0;
	while(!s.eof()) {
		std::string ln;
		std::getline(s,ln); // consume line
		if (s.eof()) break;
		std::stringstream line(ln);
		char t;
		line >> t;
		if (t=='S' || t=='s') {
			int sind;
			double sval;
			line >> sind;
			sind--;
			if (smap!=nullptr) {
				auto si = smap->find(sind);
				if (si==smap->end()) continue;
				sind = si->second;
			}
			if (!line.eof()) {
				line >> sval;
				if (sind<nsvar) tr.sx[sind] = sval;
				if (sind==wtvarnum) wt = sval;
			}
		} else if (t=='D' || t=='d') {
			int dind,dsubind;
			double dtime,dval;
			line >> dind >> dsubind >> dtime;
			dind--;
			if (dmap!=nullptr) {
				auto di = dmap->find(dind);
				if (di==dmap->end()) continue;
				dind = di->second;
			}
			if (!line.eof()) {
				dtime += (dind+dsubind/100)/100;
				dtime /= 60*60; // make units hours
				line >> dval;
				if (dind<ndvar) {
					if (tonormbywt.find(dind)!=tonormbywt.end()) {
						assert(wt>0.0);
						dval /= wt;
					}
					tr[dind].insert(dtime,dval);
					if (dtime>maxtime) maxtime = dtime;
				}
			}
		}
	}
	for(auto &vtr : tr) {
		vtr.starttime() = 0.0;
		vtr.endtime() = maxtime;
	}
}

class normmethod {
public:
	virtual ~normmethod() { }
	// age is in *months*!
	virtual double normalize(double val, double age, bool ismale) const {
		return val;
	}
	virtual double unnormalize(double val, double age, bool ismale) const {
		return val;
	}
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
	}
};
BOOST_CLASS_EXPORT_KEY(normmethod)

class linearnorm : public normmethod {
public:
	linearnorm(double mu=0.0, double std=1.0) : m(mu), s(std) {}
	linearnorm(const std::vector<traj> &ds, int varid, bool isdyn,
			int agevarnum, int sexvarnum) {
		m = 0.0; s = 1.0;
		fromdata(ds,varid,isdyn,agevarnum,sexvarnum);
	}
		
	virtual ~linearnorm() { }
	virtual double normalize(double val, double age, bool ismale) const {
		return (val-m)/s;
	}
	virtual double unnormalize(double val, double age, bool ismale) const {
		return val*s+m;
	}

protected:
	void fromdata(const std::vector<traj> &ds, int varid, bool isdyn,
			int agevarnum, int sexvarnum) {
		int n=0;
		double sum=0.0,sum2=0.0;
		if (isdyn) {
			for(auto &tr : ds) {
				if (varid>=tr.size()) continue;
				for(auto &pt : tr[varid]) { double v = normalize(pt.second.v,
						tr.sx[agevarnum],tr.sx[sexvarnum]);
					sum += v;
					sum2 += v*v;
					n++;
				}
			}
		} else {
			for(auto &tr : ds) {
				if (varid>=tr.sx.size()) continue;
				double v = this->normalize(tr.sx[varid],
						tr.sx[agevarnum],tr.sx[sexvarnum]);
				if (!std::isnan(v)) {
					sum += v;
					sum2 += v*v;
					n++;
				}
			}
		}
		if (n>0) {
			m = sum/n;
			s = sum2/n - m*m;
			if (s<=1e-10) s = 1.0;
			else s = std::sqrt(s);
		}
	}

private:
	double m,s;
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(normmethod);
		ar & BOOST_SERIALIZATION_NVP(m) & BOOST_SERIALIZATION_NVP(s);
	}
};
BOOST_CLASS_EXPORT_KEY(linearnorm)

class lookupnorm : public linearnorm {
public:
	// maps go from min age (in months) to median value
	template<typename T1, typename T2>
	lookupnorm(T1 &&femalemap, T2 &&malemap, double mu, double std)
		: linearnorm(mu,std),
			fmap(std::forward<T1>(femalemap)),
			mmap(std::forward<T2>(malemap)) {
	}
	template<typename T1, typename T2>
	lookupnorm(T1 &&femalemap, T2 &&malemap, 
			const std::vector<traj> &ds, int varid, bool isdyn,
			int agevarnum, int sexvarnum) 
		: linearnorm(0.0,1.0),
			fmap(std::forward<T1>(femalemap)),
			mmap(std::forward<T2>(malemap)) {
		fromdata(ds,varid,isdyn,agevarnum,sexvarnum);
	}
	lookupnorm() : linearnorm(0.0,1.0),
		fmap{{10000.0,1.0}}, mmap{{10000.0,1.0}} {
	}
	virtual ~lookupnorm() {}
	virtual double normalize(double val, double age, bool ismale) const {
		auto x = ismale ? mmap.upper_bound(age) : fmap.upper_bound(age);
		auto b = ismale ? mmap.begin() : fmap.begin();
		if (x!=b) --x;
		return linearnorm::normalize(val/x->second,age,ismale);
	}
	virtual double unnormalize(double val, double age, bool ismale) const {
		auto x = ismale ? mmap.upper_bound(age) : fmap.upper_bound(age);
		auto b = ismale ? mmap.begin() : fmap.begin();
		if (x!=b) --x;
		return linearnorm::unnormalize(val,age,ismale)*x->second;
	}
private:
	std::map<double,double> fmap,mmap;
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(linearnorm);
		ar & BOOST_SERIALIZATION_NVP(fmap) & BOOST_SERIALIZATION_NVP(mmap);
	}
};
BOOST_CLASS_EXPORT_KEY(lookupnorm)

typedef std::vector<shptr<normmethod>> normmap;

struct datainfo {
	bool binary;
	normmap snorm,dnorm;
	std::vector<std::string> svarnames, dvarnames;
	std::map<std::string,int> svarid,dvarid;
	std::vector<bool> sisdiscrete;
	std::vector<bool> disdiscrete;
	std::vector<std::string> episodes;
	//std::map<std::string, std::map<std::string, size_t>> sdiscretecodes;
	//std::map<std::string, std::map<std::string, size_t>> ddiscretecodes;
	//std::map<std::string, std::map<size_t, std::string>> sdiscretevals;
	//std::map<std::string, std::map<size_t, std::string>> ddiscretevals;
	std::map<std::string, std::map<std::string, size_t>> disccodes;
	std::map<std::string, std::map<size_t, std::string>> discvals;
	std::map<std::string, std::map<std::string, size_t>> contcodes;
	std::map<std::string, std::map<size_t, std::string>> contvals;
	inline std::string dvarname(int v) const {return v<0 ? "X" : dvarnames[v];}
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void save(Ar &ar, const unsigned int ver) const {
		ar & BOOST_SERIALIZATION_NVP(episodes);
		ar & BOOST_SERIALIZATION_NVP(svarnames);
		ar & BOOST_SERIALIZATION_NVP(sisdiscrete);
		//ar & BOOST_SERIALIZATION_NVP(sdiscretecodes);
		ar & BOOST_SERIALIZATION_NVP(snorm);
		ar & BOOST_SERIALIZATION_NVP(dvarnames);
		ar & BOOST_SERIALIZATION_NVP(disdiscrete);
		//ar & BOOST_SERIALIZATION_NVP(ddiscretecodes);
		ar & BOOST_SERIALIZATION_NVP(dnorm);
		ar & BOOST_SERIALIZATION_NVP(binary);
		ar & BOOST_SERIALIZATION_NVP(disccodes);
		ar & BOOST_SERIALIZATION_NVP(contcodes);
	}
	template<typename Ar>
	void load(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(episodes);
		ar & BOOST_SERIALIZATION_NVP(svarnames);
		ar & BOOST_SERIALIZATION_NVP(sisdiscrete);
		//ar & BOOST_SERIALIZATION_NVP(sdiscretecodes);
		ar & BOOST_SERIALIZATION_NVP(snorm);
		ar & BOOST_SERIALIZATION_NVP(dvarnames);
		ar & BOOST_SERIALIZATION_NVP(disdiscrete);
		//ar & BOOST_SERIALIZATION_NVP(ddiscretecodes);
		ar & BOOST_SERIALIZATION_NVP(dnorm);
		ar & BOOST_SERIALIZATION_NVP(binary);
		ar & BOOST_SERIALIZATION_NVP(disccodes);
		ar & BOOST_SERIALIZATION_NVP(contcodes);
		svarid.clear();
		for(int i=0;i<svarnames.size();i++) svarid[svarnames[i]] = i;
		dvarid.clear();
		for(int i=0;i<dvarnames.size();i++) dvarid[dvarnames[i]] = i;
/*
		sdiscretevals.clear();
		for(auto & varvalcodes : sdiscretecodes) {
			std::map<size_t, std::string> codemap;
			for(auto & valcode : varvalcodes.second)
				codemap[valcode.second] = valcode.first;
			sdiscretevals[varvalcodes.first] = codemap;
		}
		ddiscretevals.clear();
		for(auto & varvalcodes : ddiscretecodes) {
			std::map<size_t, std::string> codemap;
			for(auto & valcode : varvalcodes.second)
				codemap[valcode.second] = valcode.first;
			ddiscretevals[varvalcodes.first] = codemap;
		}
*/
		discvals.clear();
		for(auto & varvalcodes : disccodes) {
			std::map<size_t, std::string> codemap;
			for(auto & valcode : varvalcodes.second)
				codemap[valcode.second] = valcode.first;
			discvals[varvalcodes.first] = codemap;
		}
		contvals.clear();
		for(auto & varvalcodes : contcodes) {
			std::map<size_t, std::string> codemap;
			for(auto & valcode : varvalcodes.second)
				codemap[valcode.second] = valcode.first;
			contvals[varvalcodes.first] = codemap;
		}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};
BOOST_CLASS_EXPORT_KEY(datainfo)

struct dataset {
	datainfo info;
	std::vector<traj> ds;

	dataset() : info(), ds() {}
	//dataset(const dataset & rhs) {info = rhs.info; ds = rhs.ds;}

	void applynorms(int agevar, int sexvar)
	{for(traj & tr : ds) normalize(tr,agevar,sexvar);}

	void removenorms(int agevar, int sexvar)
	{for(traj &tr : ds) unnormalize(tr,agevar,sexvar);}

	void normalize(traj &tr, int agevar, int sexvar) const {
		double age = tr.sx[agevar];
		bool ismale = tr.sx[sexvar]>0.5;
		for(int i=0;i<tr.sx.size();i++) {
			tr.sx[i] = info.snorm[i]->normalize(tr.sx[i],age,ismale);
		}
		for(int i=0;i<tr.size();i++) {
			auto norm = info.dnorm[i];
			for(auto &pt : tr[i]) {
				pt.second.v = norm->normalize(pt.second.v,age,ismale);
			}
		}
	}
	void unnormalize(traj &tr, int agevar, int sexvar) const {
		double age = tr.sx[agevar];
		bool ismale = tr.sx[sexvar]>0.5;
		for(int i=0;i<tr.sx.size();i++) {
			tr.sx[i] = info.snorm[i]->unnormalize(tr.sx[i],age,ismale);
		}
		for(int i=0;i<tr.size();i++) {
			auto norm = info.dnorm[i];
			for(auto &pt : tr[i])
				pt.second.v = norm->unnormalize(pt.second.v,age,ismale);
		}
	}
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(info);
		ar & BOOST_SERIALIZATION_NVP(ds);
	}
};
BOOST_CLASS_EXPORT_KEY(dataset)

inline void akinormalize(dataset &ds) {
	auto agevari = ds.info.svarid.find("startage");
	assert(agevari!=ds.info.svarid.end());
	int agevarnum = agevari->second;
	auto sexvari = ds.info.svarid.find("sex");
	assert(sexvari!=ds.info.svarid.end());
	int sexvarnum = sexvari->second;
	
	for(unsigned int i=0;i<ds.info.svarnames.size();i++) {
		if (ds.info.snorm.size()>i) continue;
		if (ds.info.sisdiscrete[i]
				|| ds.info.svarnames[i] == "starttod"
				|| ds.info.svarnames[i] == "startage"
				|| ds.info.svarnames[i] == "baselinecreatinine")
			ds.info.snorm.emplace_back(new normmethod());
		else
			ds.info.snorm.emplace_back(new linearnorm(ds.ds,i,false,
								agevarnum,sexvarnum));
	}
	for(unsigned int i=0;i<ds.info.dvarnames.size();i++) {
		if (ds.info.dnorm.size()>i) continue;
		if (ds.info.dvarnames[i]=="Systolic Blood Pressure") {
			ds.info.dnorm.emplace_back(new lookupnorm(
				std::map<double,double>(
				{{0.0*12,86.0},{3.0*12,88.0},{4.0*12,89.0},
				 {5.0*12,91.0},{6.0*12,93.0},{7.0*12,94.0},
				 {8.0*12,96.0},{9.0*12,98.0},{10.0*12,100.0},
				 {11.0*12,102.0},{12.0*12,103.0},{13.0*12,105.0},
				 {14.0*12,107.0},{15.0*12,109.0},{16.0*12,110.0},
				 {17.0*12,111.0}}),
				std::map<double,double>(
				{{0.0*12,85.0},{3.0*12,88.0},{4.0*12,91.0},
				 {5.0*12,93.0},{6.0*12,95.0},{7.0*12,96.0},
				 {8.0*12,97.0},{9.0*12,99.0},{10.0*12,100.0},
				 {11.0*12,102.0},{12.0*12,104.0},{13.0*12,106.0},
				 {14.0*12,108.0},{15.0*12,111.0},{16.0*12,113.0},
				 {17.0*12,116.0},{18.0*12,118}}),
				ds.ds,i,true,agevarnum,sexvarnum));
		} else if (ds.info.dvarnames[i]=="Diastolic Blood Pressure") {
			ds.info.dnorm.emplace_back(new lookupnorm(
				std::map<double,double>(
				{{0.0*12,40.0},{3.0*12,45.0},{4.0*12,49.0},
				 {5.0*12,52.0},{6.0*12,54.0},{7.0*12,56.0},
				 {8.0*12,57.0},{9.0*12,58.0},{10.0*12,59.0},
				 {11.0*12,60.0},{12.0*12,61.0},{13.0*12,62.0},
				 {14.0*12,63.0},{15.0*12,64.0},{16.0*12,65.0},
				 {17.0*12,66.0}}),
				std::map<double,double>(
				{{0.0*12,37.0},{3.0*12,42.0},{4.0*12,46.0},
				 {5.0*12,50.0},{6.0*12,53.0},{7.0*12,55.0},
				 {8.0*12,57.0},{9.0*12,59.0},{10.0*12,60.0},
				 {11.0*12,61.0},{12.0*12,61.0},{13.0*12,62.0},
				 {14.0*12,62.0},{15.0*12,63.0},{16.0*12,64.0},
				 {17.0*12,65.0},{18.0*12,67}}),
				ds.ds,i,true,agevarnum,sexvarnum));
		} else if (ds.info.dvarnames[i]=="Heart Rate") {
			ds.info.dnorm.emplace_back(new lookupnorm(
				std::map<double,double>(
				{{-5.0,127.0},{0.0,143.0},{3.0,140.0},
				 {6.0,134.0},{9.0,128.0},{12.0,123.0},
				 {18.0,116.0},{24.0,110.0},{36.0,104.0},
				 {48.0,98.0},{72.0,91.0},{96.0,84.0},
				 {144.0,73.0}}),
				std::map<double,double>(
				{{-5.0,127.0},{0.0,143.0},{3.0,140.0},
				 {6.0,134.0},{9.0,128.0},{12.0,123.0},
				 {18.0,116.0},{24.0,110.0},{36.0,104.0},
				 {48.0,98.0},{72.0,91.0},{96.0,84.0},
				 {144.0,73.0}}),
				ds.ds,i,true,agevarnum,sexvarnum));
		} else if (ds.info.dvarnames[i]=="Respiratory Rate") {
			ds.info.dnorm.emplace_back(new lookupnorm(
				std::map<double,double>(
				{{0.0,43.0},{3.0,41.0},
				 {6.0,39.0},{9.0,37.0},{12.0,35.0},
				 {18.0,31.0},{24.0,28.0},{36.0,25.0},
				 {48.0,23.0},{72.0,21.0},{96.0,19.0},
				 {144.0,18.0},{180.0,16}}),
				std::map<double,double>(
				{{0.0,43.0},{3.0,41.0},
				 {6.0,39.0},{9.0,37.0},{12.0,35.0},
				 {18.0,31.0},{24.0,28.0},{36.0,25.0},
				 {48.0,23.0},{72.0,21.0},{96.0,19.0},
				 {144.0,18.0},{180.0,16}}),
				ds.ds,i,true,agevarnum,sexvarnum));
		} else ds.info.dnorm.emplace_back(new linearnorm(ds.ds,i,true,agevarnum,sexvarnum));
	}
	ds.applynorms(agevarnum,sexvarnum);
}

inline void finishload(const std::string &dir, std::ifstream &master,
		dataset &ds, int n, int skip, int svid,
		const std::map<int,int> *smap=nullptr,
		const std::map<int,int> *dmap=nullptr) {
	std::ifstream stimes((dir+"/starttimes.txt").c_str());
	int stvarid = -1;
	std::map<int,double> stimemap;
	if (stimes.good()) {
		if (svid>=0) {
			ds.info.svarnames.resize(svid+1);
			ds.info.sisdiscrete.resize(svid+1);
			ds.info.svarnames[svid] = "starttod";
			ds.info.svarid[std::string("starttod")] = svid;
			ds.info.sisdiscrete[svid] = false;
			stvarid = svid;
			svid++;
		} else {
			auto vloc = ds.info.svarid.find(std::string("starttod"));
			if (vloc!=ds.info.svarid.end()) stvarid = vloc->second;
		}
		if (stvarid!=-1) {
			while(!stimes.eof()) {
				int ind,st;
				stimes >> ind;
				if (stimes.eof()) break;
				stimes >> st;
				stimemap[ind] = (double)st/(60.0*60.0); // make hours
			}
			stimes.close();
		}
	}
	std::string ln;
	std::ifstream eplist(dir+"/eplist.txt");
	std::getline(eplist,ln); // read random seed
	std::getline(eplist,ln); // read header
	std::map<std::string,int> fn2epid;
	while(!eplist.eof()) {
		std::getline(eplist,ln);
		if (eplist.eof()) break;
		std::vector<std::string> m = strsplit(ln,',');
		fn2epid[m[1]] = std::stoi(m[0]);
	}

	std::ifstream dialist(dir+"/diagnoses/AKI HR Diagnosis Table.csv");
	std::getline(dialist,ln);
	std::vector<std::string> toks = strsplit(ln,',');
	std::vector<int> diagnum;
	for(int i=1;i<toks.size();i++) {
		std::string dname = std::string("Diag:")+toks[i];
		if (svid>=0) {
			ds.info.svarnames.resize(svid+1);
			ds.info.sisdiscrete.resize(svid+1);
			ds.info.svarnames[svid] = dname;
			ds.info.svarid[dname] = svid;
			ds.info.sisdiscrete[svid] = true;
			diagnum.emplace_back(svid);
			svid++;
		} else {
			auto vloc = ds.info.svarid.find(dname);
			if (vloc!=ds.info.svarid.end())
				diagnum.emplace_back(vloc->second);
			else diagnum.emplace_back(-1);
		}
	}
	std::map<int,std::vector<int>> diamap;
	while(!dialist.eof()) {
		std::getline(dialist,ln);
		if (dialist.eof()) break;
		std::vector<std::string> tok = strsplit(ln,',');
		std::vector<int> vals(tok.size()-1);
		int id = std::stoi(tok[0]);
		for(int i=1;i<tok.size();i++) vals[i-1] = std::stoi(tok[i]);
		diamap[id] = vals;
	}

	auto wtvari = ds.info.svarid.find("weight");
	assert(wtvari!=ds.info.svarid.end());
	int wtvarnum = wtvari->second;

	static std::set<std::string> tonormbywtstr{
		"Acyclovir","Amikacin","Amphotericin","Captopril","Cefotaxime",
		"Ceftazidime","Cidofovir","Cyclosporine","Dapsone","Enalapril",
		"Foscarnet","Ganciclovir","Gentamicin","Ibuprofen","Ketorolac",
		"Mesalamine","Methotrexate","Naproxen","","Piperacillin",
		"Sulfasalazine","Tacrolimus","Ticarcillin","Trimethoprim",
		"Valganciclovir","Vancomycin","Bumetanide","Chlorotiazide",
		"Spironolactone"};

	std::set<int> tonormbywt;
	for(auto &s : tonormbywtstr) {
		auto ll = ds.info.dvarid.find(s);
		if (ll!=ds.info.dvarid.end()) tonormbywt.insert(ll->second);
	}
		

	int id = 0;
	while(n>0 && !master.eof()) {
		std::getline(master,ln);
		if (master.eof()) break;
		if (skip>0) skip--;
		else {
			ds.ds.emplace_back(0);
			akiloadep(dir+"/"+ln,ds.ds.back(),ds.info.svarnames.size(),ds.info.dvarnames.size(),wtvarnum,tonormbywt,smap,dmap);
			if (stvarid>=0) ds.ds.back().sx[stvarid] = stimemap[id];
			auto l = fn2epid.find(ln);
			std::map<int,std::vector<int>>::const_iterator l2;
			if (l != fn2epid.end()
			    && ((l2 = diamap.find(l->second)) != diamap.end())) {
				for(int i=0;i<diagnum.size();i++)
					if (diagnum[i]>=0)
						ds.ds.back().sx[diagnum[i]] = l2->second[i];
			} else { // if not mentioned, no relevant diagnoses
				for(int i=0;i<diagnum.size();i++)
					if (diagnum[i]>=0)
						ds.ds.back().sx[diagnum[i]] = 0;
			}
			n--;
		}
		id++;
	}
	akinormalize(ds);
}

inline void akiload(std::string dir, dataset &ds, int n, int skip=0) {
	const std::set<std::string> tousesvars
		{"startage","weight","sex","admitcategory","race","baselinecreatinine"};
	std::ifstream master((dir+"/master.txt").c_str());
	if (!master.good()) return; // error out?
	std::string ln;
	bool isstatic = true;
	std::map<int,int> smap,dmap;
	int svid = 0, dvid = 0;
	while(!master.eof()) {
		std::getline(master,ln);
		if (master.eof()) return;
		if (ln == "Episodes:") break;
		if (ln == "Dynamic Variables:") isstatic = false;
		if (ln[0]=='\t') {
			if (ln[1]!='\t') {
				auto i = ln.find(": ");
				auto numstr = ln.substr(0,i);
				std::stringstream ss(ln.substr(0,i));
				int varn;
				ss >> varn;
				varn--;
				auto varname = ln.substr(i+2);
				auto pi = varname.find("(");
				if (pi!=std::string::npos) {
					while (pi>0 && varname[pi-1] == ' ') pi--;
					varname = varname.substr(0,pi);
				}
				
				if (isstatic) {
					if (tousesvars.find(varname)!=tousesvars.end()) {
						if (ds.info.svarnames.size()<=svid) {
							ds.info.svarnames.resize(svid+1);
							ds.info.sisdiscrete.resize(svid+1);
						}
						ds.info.svarnames[svid] = varname;
						ds.info.svarid[varname] = svid;
						ds.info.sisdiscrete[svid] = false;
						smap[varn] = svid;
						svid++;
					}
				} else {
					//std::cout << "dmap " << varn << " => " << varname << std::endl;
					if (ds.info.dvarnames.size()<=varn) ds.info.dvarnames.resize(varn+1);
					ds.info.dvarnames[varn] = varname;
					ds.info.dvarid[varname] = varn;
				}
			} else if (isstatic) ds.info.sisdiscrete.back() = true;
		}
	}
	finishload(dir,master,ds,n,skip,svid,&smap,nullptr);
}

inline int akicount(std::string dir) {
	std::ifstream eplist((dir+"/eplist.txt").c_str());
	if (!eplist.good()) return 0; // error out?
	std::string ln;
	std::getline(eplist,ln);
	std::getline(eplist,ln);
	int n = 0;
	while(!eplist.eof()) {
		std::getline(eplist,ln);
		if (eplist.eof()) break;
		n++;
	}
	return n;
}

inline void akiloadgiveninfo(std::string dir, dataset &ds, int n, int skip=0) {
	std::ifstream master((dir+"/master.txt").c_str());
	if (!master.good()) return; // error out?
	std::string ln;
	bool isstatic = true;
	std::map<int,int> smap,dmap;
	while(!master.eof()) {
		std::getline(master,ln);
		if (master.eof()) return;
		if (ln == "Episodes:") break;
		if (ln == "Dynamic Variables:") isstatic = false;
		if (ln[0]=='\t') {
			if (ln[1]!='\t') {
				auto i = ln.find(": ");
				auto numstr = ln.substr(0,i);
				std::stringstream ss(ln.substr(0,i));
				int varn;
				ss >> varn;
				varn--;
				auto varname = ln.substr(i+2);
				auto pi = varname.find("(");
				if (pi!=std::string::npos) {
					while (pi>0 && varname[pi-1] == ' ') pi--;
					varname = varname.substr(0,pi);
				}
				
				if (isstatic) {
					auto sloc = ds.info.svarid.find(varname);
					if (sloc!=ds.info.svarid.end())
						smap[varn] = sloc->second;
				} else {
					auto dloc = ds.info.dvarid.find(varname);
					if (dloc!=ds.info.dvarid.end())
						dmap[varn] = dloc->second;
				}
			}
		}
	}
	finishload(dir,master,ds,n,skip,-1,&smap,&dmap);
}

#endif
