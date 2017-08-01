#include "pcim.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <future>
#include <queue>
#include "serial.h"

using namespace std;

inline vector<vartrajrange> torange(const vector<traj> &data) {
	vector<vartrajrange> ret;
	if (data.empty()) return ret;
	int nv = data[0].size();
	for(auto &x : data) for(int v=0;v<nv;v++)
		ret.emplace_back(&x,v);
	return ret;
}

pcim::pcim(const vector<traj> &data,
		const vector<shptr<pcimtest>> &tests,
		const pcimparams &params) {
	const vector<vartrajrange> &d = torange(data);
	ss s = suffstats(d);
	build(d,s,tests,score(s,params),params);
}

inline vector<vartrajrange> torange(const vector<traj *> & data) {
	vector<vartrajrange> ret;
	if (data.empty()) return ret;
	int nv = data[0]->size();
	for(auto x : data) for(int v = 0; v < nv; v++) ret.emplace_back(x, v);
	return ret;
}

pcim::pcim(const vector<traj *> & data,
		const vector<shptr<pcimtest>> & tests,
		const pcimparams & params) {
	const vector<vartrajrange> & d = torange(data);
	ss s = suffstats(d);
	build(d , s, tests, score(s, params), params);
}

pcim::ss pcim::suffstats(const std::vector<vartrajrange> &data) {
	ss ret;
	ret.n=0.0;
	ret.t=0.0;
	for(const auto &x : data) {
		ret.t += x.range.second-x.range.first;
		const vartraj &vtr = (*(x.tr))[x.var];
		auto i0 = vtr.upper_bound(x.range.first);
		auto i1 = vtr.upper_bound(x.range.second);
		ret.n += distance(i0,i1);
	}
	return ret;
}

double pcim::score(const ss &d, const pcimparams &p) {
	double a_n = p.a+d.n;
	double b_n = p.b+d.t;

	double hn = d.n/2.0;
	return p.lk
		+ lgamma(a_n) - p.lga
		+ p.alb - a_n*log(b_n);
}

void pcim::calcleaf(const ss &d, const pcimparams &p) {
	rate = (p.a+d.n)/(p.b+d.t);
}

pcim::pcim(const vector<vartrajrange> &data, const pcim::ss &s,
		const vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params) {
	build(data,s,tests,basescore,params);
}

pcim::testpick pcim::picktest(const vector<vartrajrange> &data,
               const vector<shptr<pcimtest>> &tests,
			const pcimparams &params,
               int procnum, int nproc) const {
	double sc = -numeric_limits<double>::infinity();
	testpick ret;
	ret.testnum = -1; ret.s1 = ret.s2 = sc;
	for(int i=procnum;i<tests.size();i+=nproc) {
		auto &t = tests[i];
		vector<vartrajrange> td1,td2;
		for(auto &x : data) t->chop(x,td1,td2);
		if (td1.empty() || td2.empty()) continue;
		ss tss1 = suffstats(td1), tss2 = suffstats(td2);
		if (tss1.n<params.mne || tss2.n<params.mne) continue;
		double rs1 = score(tss1,params);
		double rs2 = score(tss2,params);
		if (rs1+rs2>sc) {
			ret.testnum = i;
			ret.s1 = rs1; ret.s2=rs2; sc = rs1+rs2;
			ret.ss1 = move(tss1); ret.ss2 = move(tss2);
			ret.test = t;
			ret.d1 = move(td1);
			ret.d2 = move(td2);
		}
	}
	return ret;
}

void pcim::build(const vector<vartrajrange> &data, const ss &s,
		const vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params) {
	assert(!data.empty());
	vector<vartrajrange> d1,d2;
	test.reset();
	double sc = basescore;
	testpick pick;
	if (params.nproc<=1) {
		pick = picktest(data,tests,params);
		if (pick.s1+pick.s2>sc)
			sc = pick.s1+pick.s2;
	} else {
		vector<future<testpick>> futs(params.nproc);
		for(int i=0;i<params.nproc;i++)
			futs[i] = async(launch::async,&pcim::picktest,
							this,data,tests,params,i,params.nproc);
		int picki = tests.size();
		for(auto &f : futs) {
			testpick p = f.get();
			double newsc = p.s1+p.s2;
			if (newsc>sc || (newsc==sc && p.testnum<picki)) {
				pick = move(p);
				picki = p.testnum;
				sc = newsc;
			}
		}
	}
	if (sc > basescore) {
		test = pick.test;
		ttree = shptr<pcim>(new pcim(pick.d1,pick.ss1,
								tests,pick.s1,params));
		ftree = shptr<pcim>(new pcim(pick.d2,pick.ss2,
								tests,pick.s2,params));
	} else {
		ttree.reset();
		ftree.reset();
	}
	calcleaf(s,params);
	stats = s;
}

double pcim::getevent(const traj &tr, double &t, double expsamp,
			double unisamp, int &var, double maxt) const {
	double until;
	vector<const pcim *> leaves;
	double r = getrate(tr,t,until,leaves);
	while(expsamp>(until-t)*r) {
		expsamp -= (until-t)*r;
		if (until>maxt) return maxt;
		t = until;
		r = getrate(tr,t,until,leaves);
	}
	var = leaves.size()-1;
	for(int i=0;i<leaves.size();i++) {
		unisamp -= leaves[i]->rate/r;
		if (unisamp<=0) { var = i; break; }
	}
	return t+expsamp/r;
}

double pcim::getrate(const traj &tr, double t, double &until,
			vector<const pcim *> &ret) const {
	until = numeric_limits<double>::infinity();
	ret.resize(tr.size());
	double r = 0.0;
	for(int i=0;i<tr.size();i++)
		r += getratevar(tr,i,t,until,ret[i]);
	return r;
}

double pcim::getratevar(const traj &tr, int var, double t, double &until,
			const pcim *&leaf) const {
	if (!test) { leaf = this; return rate; }
	double til;
	bool dir = test->eval(tr,var,t,til);
	if (til<until) until = til;
	return (dir ? ttree : ftree)->getratevar(tr,var,t,until,leaf);
}

void pcim::print(ostream &os) const {
	printhelp(os,0);
}
void pcim::print(ostream &os, const datainfo &info) const {
	printhelp(os,0,&info);
}

void pcim::todot(ostream &os, const datainfo &info) const {
	os << "digraph {" << endl;
	int nn = 0;
	todothelp(os,-1,false,nn,info);
	os << "}" << endl;
}

void pcim::todothelp(ostream &os, int par, bool istrue, int &nn, const datainfo &info) const {
	int mynode = nn++;
	os << "\tNODE" << mynode << " [label=\"";
	if (!ttree) os << "rate = " << rate;
	else test->print(os,info);
	os << "\"];" << endl;
	if (par>=0) os << "\tNODE" << par << " -> NODE" << mynode << " [label=\"" <<
			(istrue ? "Y: " : "N: ") << stats.n << ',' << stats.t << "\"];" << endl;
	if (ttree) {
		ttree->todothelp(os,mynode,true,nn,info);
		ftree->todothelp(os,mynode,false,nn,info);
	}
}

void pcim::printhelp(ostream &os, int lvl, const datainfo *info) const {
	for(int i=0;i<lvl;i++) os << "   ";
	if (!ttree)
		os << "rate = " << rate << " [" << stats.n << ',' << stats.t << "]" << endl;
	else {
		os << "if ";
		if (info!=nullptr) test->print(os,*info);
		else test->print(os);
		os << " [" << stats.n << ',' << stats.t << "]" << endl;
		ttree->printhelp(os,lvl+1,info);
		for(int i=0;i<lvl;i++) os << "   ";
		os << "else" << endl;
		ftree->printhelp(os,lvl+1,info);
	}
}

double pcim::llh(const std::vector<vartrajrange> &tr) const {
	if (!ttree) {
		ss d = suffstats(tr);
		return d.n*std::log(rate) - rate*d.t;
	} else {
		vector<vartrajrange> ttr,ftr;
		for(auto &x : tr) test->chop(x,ttr,ftr);
		return ttree->llh(ttr) + ftree->llh(ftr);
	}
}

BOOST_CLASS_EXPORT_IMPLEMENT(pcimtest)
BOOST_CLASS_EXPORT_IMPLEMENT(timetest)
BOOST_CLASS_EXPORT_IMPLEMENT(counttest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<counttest>)
BOOST_CLASS_EXPORT_IMPLEMENT(vartest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticgreqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticeqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(pcim)
