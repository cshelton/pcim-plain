#ifndef PCIM_H
#define PCIM_H

#include "traj.h"
#include <memory>
#include <iostream>
#include <utility>
#include <cmath>
#include <random>
#include <string>
#include "serial.h"
#include "datainfo.h"
#include <vector>
#include <array>
#include <future>

namespace boost { namespace serialization { 
	class access;
}}

typedef std::pair<double,double> timerange;

/* vartrajrange records a time range for one variable of a trajectory */
struct vartrajrange {
	vartrajrange(const traj *traject, int v, double t0, double t1) : range(t0,t1), tr(traject) {
		tr = traject;
		var = v;
	}
	vartrajrange(const vartrajrange &vtr, double t0, double t1) : range(t0,t1) {
		tr = vtr.tr;
		var = vtr.var;
	}
	vartrajrange(const traj *traject, int v)
				: range((*traject)[v].starttime(),(*traject)[v].endtime()) {
		tr = traject;
		var = v;
	}
	const traj *tr;
	int var;
	timerange range;
};

/* pcimtest is the super class of all tests (possible decisions in the
 * decision tree)
 */
class pcimtest {
public:
	virtual ~pcimtest() {} ;
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const = 0;
	// chop takes a range and breaks it according to the test
	// adds to (does not replace) outtrue and outfalse
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const = 0;
	// eval takes a trajectory, a variable (label), and a time and
	// reports the value of the test at that time for that label
	virtual bool eval(const traj &tr, int var, double t) const {
		double toss; return eval(tr,var,t,toss);
	}
	// same as above, but changed until to be the next time at which
	// the value of this test (might) change, if no other events occur
	virtual bool eval(const traj &tr, int var, double t, double &until) const = 0;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
	}
};

// tests if the raw time modulus "mod" is between "tstart" and "tend"
// stadd is a static attribute to offset the raw times
// if stadd is -1, then it is ignored
// For instance, it you want to flag morning, rush hour (and your
// time units are hours), you might set tstart = 5, tend = 9, mod = 24
// If your data all start at time 0, but this isn't really midnight,
// Let stadd be the index of a static attribute that is the actual time
// (since midnight) that a trajectory starts.
class timetest : public pcimtest {
public:
	// stadd is the index of a static variable to add to the time
	// (-1 => none)
	timetest(double tstart=0, double tend=1, double mod=1, int stadd=-1) {
		t0 = tstart; t1 = tend; m = mod; sadd = stadd;
	}
	virtual ~timetest() {}
	virtual void print(std::ostream &os) const {
		os << "time in (" << t0 << ',' << t1 << ')';
		os << " (mod " << m << "h)";
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "time in (" << t0 << ',' << t1 << ')';
		os << " (mod " << m << "h)";
	}
	
	inline double remerge(double base, double inc) const {
		return base+inc;
	}

	inline double breakup(double v, double &n) const {
		double ret = std::remainder(v,m);
		if (ret<0) ret += m;
		n = v - ret;
		return ret;
	}
			
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		int igt = outtrue.size();
		int igf = outfalse.size();
		double temp;
		double del = sadd<0 ? 0 : -in.tr->sx[sadd];
		double myt0 = breakup(t0+del,temp);
		double myt1 = breakup(t1+del,temp);
		double startbase,endbase;
		double start = breakup(in.range.first,startbase);
		double end = breakup(in.range.second,endbase);
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		std::vector<vartrajrange> &out0 = myt0<myt1 ? outfalse : outtrue;
		std::vector<vartrajrange> &out1 = myt0<myt1 ? outtrue: outfalse;
		auto startpt = in.range.first; //remerge(startbase,start);
		while(startbase<endbase) {
			if (start<tmin) {
				auto nextpt = remerge(startbase,tmin);
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				startpt = nextpt;
				start = tmin;
			}
			if (start<tmax) {
				auto nextpt = remerge(startbase,tmax);
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				startpt = nextpt;
				start = tmax;
			}
			start -= m;
			startbase += m;
		}
		if (start<tmin) {
			if (end<tmin) {
				auto nextpt = in.range.second;
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				return;
			} else {
				auto nextpt = remerge(startbase,tmin);
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				start = tmin;
				startpt = nextpt;
			}
		}
		if (start<tmax) {
			if (end<tmax) {
				auto nextpt = in.range.second;
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				return;
			} else {
				auto nextpt = remerge(startbase,tmax);
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				start = tmax;
				startpt = nextpt;
			}
		}
		if (start<end) {
			auto nextpt = in.range.second;
			if (startpt<nextpt)
				out0.emplace_back(in,startpt,nextpt);
		}
	}

	virtual bool eval(const traj &tr, int var, double t) const {
		double temp;
		double del = sadd<0 ? 0 : -tr.sx[sadd];
		double myt0 = breakup(t0+del,temp);
		double myt1 = breakup(t1+del,temp);
		double tbase;
		double tmod = breakup(t,tbase);
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		if (tmod<tmin || tmod>tmax) return myt0>myt1;
		return myt1>myt0;
	}
	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		double temp;
		double del = sadd<0 ? 0 : -tr.sx[sadd];
		double myt0 = breakup(t0+del,temp);
		double myt1 = breakup(t1+del,temp);
		double tbase;
		double tmod = breakup(t,tbase);
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		if (tmod<tmin) {
			until = remerge(tbase,tmin);
			assert(until>t);
			return myt0>myt1;
		}
		if (tmod<tmax) {
			until = remerge(tbase,tmax);
			assert(until>t);
			return myt1>myt0;
		}
		until = remerge(tbase,tmin)+m;
		assert(until>t);
		return myt0>myt1;
	}

private:
	double t0,t1,m;
	int sadd;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(t0) & BOOST_SERIALIZATION_NVP(t1) & BOOST_SERIALIZATION_NVP(m) & BOOST_SERIALIZATION_NVP(sadd);
	}
};

template<typename D>
class varstattest : public pcimtest {
public:
	varstattest(int testvar=0, double lag0=0, double lag1=1) {
		v = testvar;
		maxlag = std::max(lag0,lag1);
		minlag = std::min(lag0,lag1);
	}
	virtual ~varstattest() {}
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const =0;
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const vartraj &tr = (*(in.tr))[v==-1?in.var:v];
		const auto &e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0-maxlag);
		auto i1 = tr.upper_bound(t0-minlag);
		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
		double t1 = t0;
		bool currval = static_cast<const D *>(this)->evalstat(stat);
		while(t1<tend && i0!=e) {
			t1 = std::min(i0==e ?
					std::numeric_limits<double>::infinity() : i0->first+maxlag,
					i1==e ?
					std::numeric_limits<double>::infinity() : i1->first+minlag);
			while (i0!=e && i0->first+maxlag<=t1) {
				stat.del(i0->first,i0->second);
				++i0;
			}
			while (i1!=e && i1->first+minlag<=t1) {
				stat.add(i1->first,i1->second);
				++i1;
			}
			if (t1>=tend) t1 = tend;
			bool newval = static_cast<const D *>(this)->evalstat(stat);
			if (t0<t1) {
				if (!newval && currval) {
					outtrue.emplace_back(in,t0,t1);
					t0 = t1;
					currval = false;
				} else if (newval && !currval) {
					outfalse.emplace_back(in,t0,t1);
					t0 = t1;
					currval = true;
				}
			}
		}
		if (t0<tend) {
			if (currval) outtrue.emplace_back(in,t0,tend);
			else outfalse.emplace_back(in,t0,tend);
		}

	}
	virtual bool eval(const traj &tr, int var, double t) const {
		const vartraj &vtr = tr[v==-1?var:v];
		double tnext
			= std::nextafter(tnext,std::numeric_limits<double>::infinity());
		double t0 = tnext-maxlag;
		if (t0==t-maxlag)
			t0 = std::nextafter(t0,std::numeric_limits<double>::infinity());
		double t1 = tnext-minlag;
		if (t1==t-minlag)
			t1 = std::nextafter(t1,std::numeric_limits<double>::infinity());
		auto i0 = vtr.lower_bound(t0);
		auto i1 = vtr.lower_bound(t1);

		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
		return static_cast<const D *>(this)->evalstat(stat);
	}

	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		const vartraj &vtr = tr[v==-1?var:v];
		const auto &e = vtr.cend();
		double tnext
			= std::nextafter(t,std::numeric_limits<double>::infinity());
		double t0 = tnext-maxlag;
		if (t0==t-maxlag)
			t0 = std::nextafter(t0,std::numeric_limits<double>::infinity());
		double t1 = tnext-minlag;
		if (t1==t-minlag)
			t1 = std::nextafter(t1,std::numeric_limits<double>::infinity());
		auto i0 = vtr.lower_bound(t0);
		auto i1 = vtr.lower_bound(t1);
		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
		until = std::min(i0!=e ? i0->first+maxlag
				: std::numeric_limits<double>::infinity(),
				i1!=e ? i1->first+minlag
				: std::numeric_limits<double>::infinity());
		assert(until>t);
		return static_cast<const D *>(this)->evalstat(stat);
	}

protected:
	double minlag,maxlag;
	int v;

private:
	friend class boost::serialization::access;
	template<typename Ar>
		void serialize(Ar &ar, const unsigned int ver) {
			ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
			ar & BOOST_SERIALIZATION_NVP(minlag) & BOOST_SERIALIZATION_NVP(maxlag) & BOOST_SERIALIZATION_NVP(v);
		}
};


// test if count of number of events of var testvar (-1 == currvar)
// from t-lag0 to t-lag1 is greater than thresh
class counttest : public varstattest<counttest> {
public:
	counttest(int thresh=0, int testvar=0, double lag0=0, double lag1=1)
			: varstattest<counttest>(testvar,lag0,lag1) { theta=thresh; }
	virtual ~counttest() {}
	virtual void print(std::ostream &os) const {
		os << "# " << v << " in [" << maxlag << ',' << minlag << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		//os << "# " << info.dvarname(v) << " measurements in [" 
		//	<< maxlag << ',' << minlag << ") >= " << theta;
		os << "# " << info.dvarname(v);
		os << " measrmnts in last " 
			<< maxlag << " time >= " << theta;
	}

	struct statT {
		statT() { n=0; }
		int n;
		void add(double,vartraj::rec) { n++; }
		void del(double,vartraj::rec) { n--; }
	};

	bool evalstat(const statT &s) const {
		return s.n>=theta;
	}
	
protected:
	int theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		//ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(varstattest<counttest>);
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<counttest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta);
	}
};

// test is current variable == testvar
class vartest : public pcimtest {
public:
	vartest(int testvar=0) : pcimtest() { v = testvar; };
	virtual ~vartest() {} ;
	virtual void print(std::ostream &os) const { os << "var == " << v; }
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "X == " << info.dvarname(v);
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.var==v) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const traj &tr, int var, double t) const {
		return var==v;
	}
	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		return var==v;
	}
private:
	int v;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v);
	}
};

// tests a static value against testval.
class staticgreqtest : public pcimtest {
public:
	staticgreqtest(double testval=0, int testvar=0) : pcimtest() {
		v = testvar; theta=testval;
	};
	virtual ~staticgreqtest() {} ;
	virtual void print(std::ostream &os) const {
		os << "svar(" << v << ") >= " << theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << info.svarnames[v] << " >= " << theta;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.tr->sx[v]>=theta) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const traj &tr, int var, double t) const {
		return tr.sx[v]>=theta;
	}
	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		return tr.sx[v]>=theta;
	}
private:
	int v;
	double theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v) & BOOST_SERIALIZATION_NVP(theta);
	}
};

// tests a static value against testval.
class staticeqtest : public pcimtest {
public:
	staticeqtest(double testval=0, int testvar=0) : pcimtest() {
		v = testvar; theta=testval;
	};
	virtual ~staticeqtest() {} ;
	virtual void print(std::ostream &os) const {
		os << "svar(" << v << ") == " << theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << info.svarnames[v] << " == " << theta;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.tr->sx[v]==theta) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const traj &tr, int var, double t) const {
		return tr.sx[v]==theta;
	}
	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		return tr.sx[v]>=theta;
	}
private:
	int v;
	double theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v) & BOOST_SERIALIZATION_NVP(theta);
	}
};

class pcim {
public:
	class pcimparams {
	friend pcim;
	public:
		// alpha -- pseudo event counts
		// beta -- pseudo durations
		// kappa -- kappa^n is prior on tree of size n (kappa<=1!)
		// minnumevents -- minimum number of events per leaf
		// numproc -- max number of threads to start while learning
		pcimparams(const pcimparams & p) : a(p.a), b(p.b), k(p.k),
				mne(p.mne), nproc(p.nproc) { set(); }
		pcimparams(double alpha, double beta, double kappa,
				int minnumevents=0, int numproc=1) {
			a = alpha; b = beta; k = kappa;
			nproc = numproc;
			mne = minnumevents;
			set();
		}
	private:
		double a,b,k;
		double alb; // alpha*log(beta)
		double lga; // log(gamma(alpha))
		double lk; // log(kappa)
		int nproc,mne;
		void set() {
			alb = a * ::log(b);
			lga = std::lgamma(a);
			lk = ::log(k);
		}
	};

	virtual ~pcim() {}

	// learn a pcim (c-pcim, really) from data
	// tests is a set of possible tests
	pcim(const std::vector<traj> & data,
		const std::vector<shptr<pcimtest>> & tests,
		const pcimparams & params);

	// same, but traj are ptrs
	pcim(const std::vector<traj *> & data,
		const std::vector<shptr<pcimtest>> & tests,
		const pcimparams & params);

	// build from subtrees
	pcim(shptr<pcimtest> tst, shptr<pcim> truebranch, shptr<pcim> falsebranch)
		: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0;
	}

	// same, but with ptrs
	pcim(pcimtest * tst, pcim * truebranch, pcim * falsebranch)
		: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0;
	}

	pcim(double lambda = 1.0) 
		: ttree(), ftree(), test() {
		rate =lambda;
	}

	// generates the end of a trajectory
	template<typename R>
	double samplecomplete(traj & ret, double T, R & rand) const {
		double t = 0.0;
		int var;
		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0, 1.0);
		double lastt = t;
		while((t=getevent(ret,lastt,expdist(rand),unifdist(rand),var,T))<T) {
			ret[var].insert(t);
			lastt = t;
		}
		return lastt;
	}

	template<typename R>
	traj sample(double T,int nvar,R &rand) const {
		traj ret(nvar);
		if (T<0.0) return ret;
		for(int i=0;i<nvar;i++) {
			ret[i].starttime() = 0.0;
			ret[i].endtime() = T;
		}
		samplecomplete(ret,T,rand);
		return ret;
	}

	// returns relevant leaves in ret and sum as return value
	double getrate(const traj &tr, double t, double &until, std::vector<const pcim *> &ret) const;
	// returns new time and sets var to the variable (label)
	double getevent(const traj &tr, double &t, double expsamp,
			double unisamp, int &var, double maxt) const;

	// these output the tree in a semi-readable form
	// "todot" outputs in a format readable by graphviz
	void print(std::ostream &os) const;
	void print(std::ostream &os, const datainfo &info) const;
	void todot(std::ostream &os, const datainfo &info) const;

	void save(std::ostream &os) const;
	void load(std::ostream &os);

	double llh(const traj &tr) const {
		std::vector<vartrajrange> vtr;
		for(int v=0;v<tr.size();v++) vtr.emplace_back(&tr,v);
		return llh(vtr);
	}

	inline void printtest(std::ostream & os, const datainfo & info) const {
		test->print(os, info);
	}

private:
	class ss {
	public:
		double n; // not int, in case need expected value
		double t; // total time
		ss() : n(0.0), t(0.0) {}
	private:
		friend class boost::serialization::access;
		template<typename Ar>
		void serialize(Ar &ar, const unsigned int ver) {
			ar & BOOST_SERIALIZATION_NVP(n) & BOOST_SERIALIZATION_NVP(t);
		}
	};

	pcim(const std::vector<vartrajrange> &data, const pcim::ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params);

	void build(const std::vector<vartrajrange> &data, const ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params);
	
	static ss suffstats(const std::vector<vartrajrange> &data);
	static double score(const ss &d, const pcimparams &p);
	void calcleaf(const ss &d, const pcimparams &p);

	double getratevar(const traj &tr, int var, double t, double &until,
				const pcim *&leaf) const;

	void printhelp(std::ostream &os, int lvl, const datainfo *info=nullptr) const;
	void todothelp(std::ostream &os, int par, bool istrue, int &nn, const datainfo &info) const;

	struct testpick {
		int testnum;
		shptr<pcimtest> test;
		std::vector<vartrajrange> d1,d2;
		double s1,s2;
		ss ss1,ss2;
	};

	testpick picktest(const std::vector<vartrajrange> &data,
			const std::vector<shptr<pcimtest>> &tests,
			const pcimparams &params,
			int procnum=0, int nproc=1) const;



	void trajtofeatures(const std::vector<vartrajrange> &tr,
				std::vector<double> &f) const;

	std::vector<double> trajtofeatures(
			const traj & tr,
			const std::vector<double> & endtimes,
			const size_t states) const;
	std::vector<double> varstatedurations(const vartrajrange & vtr,
		const size_t states) const;
	typedef std::pair<vartraj::const_iterator, vartraj::const_iterator> itemT;
	class eventcomp {
	public:
		eventcomp(const bool reverse = false) : reverse(reverse) {}
		bool operator ()(const itemT & a, const itemT & b) const {
			if (reverse) return a.first->first > b.first->first;
			else return a.first->first < b.first->first;
		}
	private:
		bool reverse;
	};


	double similarity(const std::vector<vartrajrange> &tr1,
				const std::vector<vartrajrange> &tr2) const;

	shptr<pcim> ftree,ttree;
	shptr<pcimtest> test;
	double rate;
	ss stats;
public:
	
private:
	double llh(const std::vector<vartrajrange> &tr) const;

	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
	}
};

BOOST_CLASS_EXPORT_KEY(pcimtest)
BOOST_CLASS_EXPORT_KEY(timetest)
BOOST_CLASS_EXPORT_KEY(counttest)
BOOST_CLASS_EXPORT_KEY(vartest)
BOOST_CLASS_EXPORT_KEY(staticgreqtest)
BOOST_CLASS_EXPORT_KEY(staticeqtest)
BOOST_CLASS_EXPORT_KEY(pcim)
BOOST_SERIALIZATION_SHARED_PTR(pcimtest)
BOOST_SERIALIZATION_SHARED_PTR(pcim)

#endif
