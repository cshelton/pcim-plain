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
#include "load.h"
#include <vector>
#include <array>
#include <future>

// Use Fisher information to compute the features
#define USEINFINV

// if defined, uses "persistance arcs" (ie next value depends on previous)
#define USEPERSIST

// Extract features using only rate
//#define RATEONLY

#ifdef USEPERSIST
#include "Eigen/Dense"
#endif

namespace boost { namespace serialization { 
	class access;
}}

typedef std::pair<double,double> timerange;

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

class pcimtest {
public:
	virtual ~pcimtest() {} ;
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const = 0;
	// adds to (does not replace) outtrue and outfalse
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const = 0;
	virtual bool eval(const traj &tr, int var, double t) const {
		double toss; return eval(tr,var,t,toss);
	}
	virtual bool eval(const traj &tr, int var, double t, double &until) const = 0;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
	}
};

// tests if last value of testvar >= thresh (if no last value, value taken
//  to be 0)
class lasttest : public pcimtest {
public:
	lasttest(int testvar=0, double thresh=0) {
		v = testvar;
		theta = thresh;
	}
	virtual ~lasttest() {}
	virtual void print(std::ostream &os) const {
		os << "most recent " << v << " >= " << theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "most recent value for " << info.dvarname(v) << " >= " << theta;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const vartraj &tr = (*(in.tr))[v==-1?in.var:v];
		const auto e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0);
		if (i0!=tr.begin()) --i0;
		auto i1 = tr.lower_bound(tend);
		bool currval = (i0==e ? 0.0 : i0->second.v) >=theta;
		double t1 = i0->first;
		while(t1<tend && i0!=e) {
			++i0;
			while (i0!=e && i0->first<=t1) ++i0;
			if (i0==e || i0->first>=tend) t1 = tend;
			else t1 = i0->first;
			bool newval = (i0==e ? 0.0 : i0->second.v) >= theta;
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
		if (t0<tend) {
			if (currval) outtrue.emplace_back(in,t0,tend);
			else outfalse.emplace_back(in,t0,tend);
		}
	}
	virtual bool eval(const traj &tr, int var, double t) const {
		const vartraj &vtr = tr[v==-1?var:v];
		if (vtr.empty()) return 0.0 >= theta;
		auto i0 = vtr.lower_bound(t);
		if (i0==vtr.cend() || i0->first>t) --i0;
		return (i0==vtr.cend() ? 0.0 : i0->second.v) >= theta;
	}

	virtual bool eval(const traj &tr, int var, double t, double &until) const {
		const vartraj &vtr = tr[v==-1?var:v];
		if (vtr.empty()) {
			until = std::numeric_limits<double>::infinity();
			return 0.0 >= theta;
		}
		auto i0 = vtr.lower_bound(t);
		auto e = vtr.cend();
		auto i1 = i0;
		if (i0==e || i0->first>t) --i0;
		else ++i1;
		until = i1!=e ? i1->first : std::numeric_limits<double>::infinity();
		assert(until>t);
		return (i0==e ? 0.0 : i0->second.v) >= theta;
	}

protected:
	double theta;
	int v;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v) & BOOST_SERIALIZATION_NVP(theta);
	}
};

class timetest : public pcimtest {
public:
	// stadd is the index of a static variable to add to the time
	// (-1 => none)
	timetest(double tstart=0, double tend=1, double mod=1, int stadd=-1) {
		t0 = tstart; t1 = tend; m = mod; sadd = stadd;
	}
	virtual ~timetest() {}
	virtual void print(std::ostream &os) const {
		//os << "time (mod " << m/3600 << "h) in (" << t0/3600 << ',' << t1/3600 << ')';
		os << "time in (" << t0/3600 << ',' << t1/3600 << ')';
		os << " (mod " << m/3600 << "h)";
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		//os << "time (mod " << m/3600 << "h) in (" << t0/3600 << ',' << t1/3600 << ')';
		os << "time in (" << t0/3600 << ',' << t1/3600 << ')';
		os << " (mod " << m/3600 << "h)";
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
				//dump(in,outtrue,outfalse,igt,igf);
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
				//dump(in,outtrue,outfalse,igt,igf);
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
		//dump(in,outtrue,outfalse,igt,igf);
	}

	void dump(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse,
			int igt, int igf) const {
		std::cout << "mod" << m << " in (" << t0 <<','<< t1 <<")" << std::endl;
		std::cout << "in: " << in.range.first << ' ' << in.range.second << std::endl;
		std::cout << "outtrue: " << std::endl;
		for(int i=igt;i<outtrue.size();i++) {
			auto &x = outtrue[i];
			std::cout << "\t" << x.range.first << ' ' << x.range.second << std::endl;
		}
		std::cout << "outfalse: " << std::endl;
		for(int i=igf;i<outfalse.size();i++) {
			auto &x = outfalse[i];
			std::cout << "\t" << x.range.first << ' ' << x.range.second << std::endl;
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
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second.v);
		double t1 = t0;
		bool currval = static_cast<const D *>(this)->evalstat(stat);
		while(t1<tend && i0!=e) {
			t1 = std::min(i0==e ?
				std::numeric_limits<double>::infinity() : i0->first+maxlag,
						i1==e ?
				std::numeric_limits<double>::infinity() : i1->first+minlag);
			while (i0!=e && i0->first+maxlag<=t1) {
				stat.del(i0->first,i0->second.v);
				++i0;
			}
			while (i1!=e && i1->first+minlag<=t1) {
				stat.add(i1->first,i1->second.v);
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
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second.v);
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
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second.v);
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
		os << "# " << v << " in [" << maxlag/3600 << ',' << minlag/3600 << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		//os << "# " << info.dvarname(v) << " measurements in [" 
		//	<< maxlag/3600 << ',' << minlag/3600 << ") >= " << theta;
		os << "# ";
		if (info.dvarname(v) == "Systolic blood pressure non-invasive (mmHg)")
			os << "Systolic BP non-inv";
		else if (info.dvarname(v) == "Pulse oximetry (%)")
			os << "Pulse ox (%)";
		else if (info.dvarname(v) == "Temperature ( C)")
			os << "Temp (C)";
		else
			os << info.dvarname(v);
		os << " measrmnts in last " 
			<< maxlag/3600 << " hours >= " << theta;
	}

	struct statT {
		statT() { n=0; }
		int n;
		void add(double,double) { n++; }
		void del(double,double) { n--; }
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

// test if mean of values of var testvar (-1 == currvar)
// from t-lag0 to t-lag1 is greater than thresh
// if no value in interval, mean is taken to be 0.0
class meantest : public varstattest<meantest> {
public:
	meantest(double thresh=0, int testvar=0, double lag0=0, double lag1=1)
			: varstattest<meantest>(testvar,lag0,lag1) { theta=thresh; }
	virtual ~meantest() {}
	virtual void print(std::ostream &os) const {
		os << "mean " << v << " in [" << maxlag/3600 << ',' << minlag/3600 << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "mean " << info.dvarname(v) << " in [" << maxlag/3600 << ',' << minlag/3600 << ") >= "
				<< theta;
	}

	struct statT {
		statT() { n=0; m=0.0; }
		int n;
		double m;
		void add(double,double v) { n++; m += v;}
		void del(double,double v) { n--; m -= v;}
	};

	bool evalstat(const statT &s) const {
		return (s.n >0 ? s.m/s.n : 0.0) >=theta;
	}
	
protected:
	double theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<meantest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta);
	}
};

// test if variance of values of var testvar (-1 == currvar)
// from t-lag0 to t-lag1 is greater than thresh
// if no value in interval, variance is taken to be 0.0
class stdvtest : public varstattest<stdvtest> {
public:
	stdvtest(double thresh=0, int testvar=0, double lag0=0, double lag1=1)
			: varstattest<stdvtest>(testvar,lag0,lag1) { theta=thresh; }
	virtual ~stdvtest() {}
	virtual void print(std::ostream &os) const {
		os << "variance " << v << " in [" << maxlag/3600 << ',' << minlag/3600 << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "variance " << info.dvarname(v) << " in [" << maxlag/3600 << ',' << minlag/3600 << ") >= "
				<< theta;
	}

	struct statT {
		statT() { n=0; m=0.0; m2 = 0.0;}
		int n;
		double m;
		double m2;
		void add(double,double v) { n++; m += v; m2 += v * v;}
		void del(double,double v) { n--; m -= v; m2 -= v * v;}
	};

	bool evalstat(const statT &s) const {
		return (s.n > 0 ? (s.m2 - s.m * s.m / s.n) / s.n : 0.0) >= theta;
	}
	
protected:
	double theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<stdvtest>>(*this));
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
#ifdef USEPERSIST
	static constexpr size_t npredfeat = 2;
	// below should go into params at some point
	static constexpr double invarreg = 100;

	typedef Eigen::Matrix<double,npredfeat,1> wtT;
	typedef Eigen::Matrix<double,npredfeat,npredfeat> wtvarT;
	static const wtT zerowt; // = wtT::Zero(npredfeat,1);
	static const wtvarT zerowtvar; // = wtT::Zero(npredfeat,npredfeat);

	inline static double quad(const wtT &x, const wtvarT &M) {
		return x.transpose()*M*x;
	}
	inline static wtT div(const wtT &x, const wtvarT &M) {
		return M.ldlt().solve(x);
	}
#else
	typedef double wtT;
	typedef double wtvarT;
	static constexpr double zerowt = 0.0;
	static constexpr double zerowtvar = 0.0;
	inline static double quad(const wtT &x, const wtvarT &M) {
		return x*M*x;
	}
	inline static wtT div(const wtT &x, const wtvarT &M) {
		return x/M;
	}
#endif
	class pcimparams {
	friend pcim;
	public:
		// alpha -- pseudo event counts
		// beta -- pseudo durations
		// kappa -- kappa^n is prior on tree of size n (kappa<=1!)
		// valalpha, valbeta, valkappa, valmu -- normal-gamma
		// 		priors for Gaussian distribution on event *value*
		pcimparams(const pcimparams & p) : a(p.a), b(p.b), k(p.k), va(p.va),
		vb(p.vb), vk(p.vk), vm(p.vm), mne(p.mne), nproc(p.nproc) {set();}
		pcimparams(double alpha, double beta, double kappa,
				const wtT &valmu, const wtvarT &valkappa, double valalpha,
				double valbeta, int minnumevents=0, int numproc=1) {
			a = alpha; b = beta; k = kappa;
			va = valalpha; vb = valbeta; vk = valkappa; vm = valmu;
			nproc = numproc;
			mne = minnumevents;
			set();
		}
	private:
		double a,b,k;
		double va,vb;
		wtvarT vk;
		wtT vm;
		double alb; // alpha*log(beta)
		double lga; // log(gamma(alpha))
		double lk; // log(kappa)
		double lgva; // log(gamma(valalpha))
		double valvb; // valalpha * log(valbeta)
		double lvk; // log(|valkappa|)
		double m2k; // valmu^2*valkappa
		int nproc,mne;
		void set() {
			alb = a * ::log(b);
			valvb = va * ::log(vb);
			lga = std::lgamma(a);
			lgva = std::lgamma(va);
			lk = ::log(k);
			m2k = pcim::quad(vm,vk);
#ifdef USEPERSIST
			lvk = ::log(vk.determinant());
#else
			lvk = ::log(vk);
#endif
		}
	};

	virtual ~pcim() {}

	pcim(
		const std::vector<traj> & data,
		const std::vector<shptr<pcimtest>> & tests,
		const pcimparams & params, size_t trackbranch = 0);

	pcim(
		const std::vector<traj *> & data,
		const std::vector<shptr<pcimtest>> & tests,
		const pcimparams & params, size_t trackbranch = 0);

	pcim(shptr<pcimtest> tst, shptr<pcim> truebranch, shptr<pcim> falsebranch)
		: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0;
#ifdef USEPERSIST
		mu = wtT::Zero(npredfeat, 1);
#else
		mu = 0.0;
#endif
		sigma = 1.0;
	}

	pcim(pcimtest * tst, pcim * truebranch, pcim * falsebranch)
		: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0;
#ifdef USEPERSIST
		mu = wtT::Zero(npredfeat, 1);
#else
		mu = 0.0;
#endif
		sigma = 1.0;
	}

	pcim(double lambda = 1.0, wtT mean = zerowt, double stddev = 1.0)
		: ttree(), ftree(), test() {
		rate =lambda;
		mu = mean;
		sigma = stddev;
		globalm = 0;
	}

	template<typename R>
	double samplecomplete(traj & ret, double T, R & rand) const {
		double t = 0.0;
		int var;
		double val;
		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0, 1.0);
		std::normal_distribution<> normdist(0.0, 1.0);
		double lastt = t;
		while((t = getevent(ret, lastt, expdist(rand), unifdist(rand),
							normdist(rand), var, val, T)) < T) {
			ret[var].insert(t, val);
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

	template<typename Tst, typename R>
	std::pair<double,double> eventprobhelp(traj tr, double T, int vnum, Tst &&tst, R &rand, int nsamp=100, int nsubsample=5000) const {
		std::normal_distribution<> normdist(0.0,1.0);
		double truewt = 0.0, wt = 0.0;
		for(int i=0;i<nsamp;i++) { // replace nsamp with bound on accuracy!
			traj h{tr};
			double lastt = samplecomplete(h,T,rand);
			double until,q;
			const pcim *leaf;
			q = getratevar(h,vnum,T,until,leaf);
			double w = q; //*exp(-q*(T-lastt));
			wt += w*nsubsample;
			for(int i=0;i<nsubsample;i++) {
#ifdef USEPERSIST
				double x = 0.0;
				if (!tr[vnum].empty()) x = tr[vnum].rbegin()->second.v;
				double y = predfeat(x).transpose()*leaf->mu
						+ normdist(rand)*leaf->sigma ;
#else
				double y = leaf->mu + normdist(rand)*leaf->sigma ;
#endif
				if (tst(y)) truewt += w;
			}
		
		}
		return std::make_pair(truewt,wt);
	}

	template<typename Tst, typename R>
	double eventprob(traj tr, double T, int vnum, Tst &&tst, R &rand, int nproc=1, int nsamp=100, int nsubsample=5000) const {
		double truewt = 0.0, wt = 0.0;
		if (nproc<=1) std::tie(truewt,wt) = eventprobhelp(tr,T,vnum,tst,rand,nsamp,nsubsample);
		else {
			std::vector<std::future<std::pair<double,double>>> futs(nproc);
			std::vector<R> rs;
			for(int i=0;i<nproc;i++) 
				rs.emplace_back(rand());
			for(int i=0;i<nproc;i++) {
				auto &rr = rs[i];
				futs[i] = async(std::launch::async,
					[this,&tr,T,vnum,tst,&rr,nsubsample](int n) {
						return this->eventprobhelp(tr,T,vnum,tst,rr,n,nsubsample);
					},
					nsamp*(i+1)/nproc - nsamp*i/nproc);
			}
			double ltruewt, lwt;
			for(auto &f : futs) {
				std::tie(ltruewt,lwt) = f.get();
				truewt += ltruewt;
				wt += lwt;
			}
		}
		return truewt/wt;
	}

	template<typename R>
	std::pair<double,double> eventprobthreshhelp(const traj &tr, double T, int vnum,
						double thresh, R &rand, int nsamp=100) const {
		std::normal_distribution<> normdist(0.0,1.0);
		double truewt = 0.0, wt = 0.0;
		const double sqrt2 = std::sqrt(2.0);
		for(int i=0;i<nsamp;i++) { // replace nsamp with bound on accuracy!
			traj h{tr};
			double lastt = samplecomplete(h,T,rand);
			double until,q;
			const pcim *leaf;
			q = getratevar(h,vnum,T,until,leaf);
			double w = q; //*exp(-q*(T-lastt));
			wt += w;
#ifdef USEPERSIST
			double x = 0.0;
			if (!tr[vnum].empty()) x = tr[vnum].rbegin()->second.v;
			double meany = predfeat(x).transpose()*leaf->mu;
#else
			double meany = leaf->mu;
#endif
			truewt += w*(0.5-std::erf((thresh-meany)/(leaf->sigma*sqrt2))/2.0);
		}
		return std::make_pair(truewt,wt);
	}

	template<typename R>
	double eventprobthresh(const traj &tr, double T, int vnum,
						double thresh, R &rand, int nproc=1, int nsamp=100) const {
		double truewt = 0.0, wt = 0.0;
		if (nproc<=1) std::tie(truewt,wt) = eventprobthreshhelp(tr,T,vnum,thresh,rand,nsamp);
		else {
			std::vector<std::future<std::pair<double,double>>> futs(nproc);
			std::vector<R> rs;
			for(int i=0;i<nproc;i++)
				rs.emplace_back(rand());
			for(int i=0;i<nproc;i++) {
				auto &rr = rs[i];
				futs[i] = async(std::launch::async,
					//eventprobthreshhelp,this,tr,T,vnum,thresh,rr,
					[this,&tr,T,vnum,thresh,&rr](int n) {
						return this->eventprobthreshhelp(tr,T,vnum,thresh,rr,n);
					},
					nsamp*(i+1)/nproc - nsamp*i/nproc);
			}
			double ltruewt, lwt;
			for(auto &f : futs) {
				std::tie(ltruewt,lwt) = f.get();
				truewt += ltruewt;
				wt += lwt;
			}
		}
		return truewt/wt;
	}

	// returns relevant leaves in ret and sum as return value
	double getrate(const traj &tr, double t, double &until, std::vector<const pcim *> &ret) const;
	// returns new time and sets var and val to the variable and its value
	double getevent(const traj &tr, double &t, double expsamp, double unisamp, double normsamp,
					int &var, double &val, double maxt) const;

	void print(std::ostream &os) const;
	void print(std::ostream &os, const datainfo &info) const;
	void todot(std::ostream &os, const datainfo &info) const;

	void save(std::ostream &os) const;
	void load(std::ostream &os);

	std::vector<std::string> featurenames() const {
		std::vector<std::string> ret;
		featurenames(ret,"");
		return ret;
	}
	std::vector<double> trajtofeatures(const traj &tr) const {
		std::vector<double> ret;
		std::vector<vartrajrange> vtr;
		for(int v=0;v<tr.size();v++) vtr.emplace_back(&tr,v);
		trajtofeatures(std::vector<vartrajrange>{vtr},ret);
		return ret;
	}
	double llh(const traj &tr) const {
		std::vector<vartrajrange> vtr;
		for(int v=0;v<tr.size();v++) vtr.emplace_back(&tr,v);
		return llh(vtr);
	}

	std::vector<double> trajtofeatures(
			const traj & tr, const size_t states) const {
		std::vector<double> endtimes;
		for(size_t v = 0; v < tr.size(); v++)
			endtimes.emplace_back(tr[v].endtime());
const double endtime = tr[0].endtime();
for (const auto & vt : tr) assert(vt.endtime() == endtime);
for (const auto & vt : tr) assert(vt.starttime() == 0.0);
		return trajtofeatures(tr, endtimes, states);
	}
	double similarity(const traj &tr1, const traj &tr2) const {
		std::vector<vartrajrange> vtr1;
		for(int v=0;v<tr1.size();v++)
			vtr1.emplace_back(&tr1,v);
		std::vector<vartrajrange> vtr2;
		for(int v=0;v<tr2.size();v++)
			vtr2.emplace_back(&tr2,v);
		return similarity(std::vector<vartrajrange>{vtr1},
					std::vector<vartrajrange>{vtr2});
	}

	void id_leaves(const shptr<pcim> & node, size_t & id);
	inline bool is_leaf() const {return test == 0;}
	inline shptr<pcim> tree(const std::string & branch) const {
		if(branch == "T") return ttree;
		if(branch == "F") return ftree;
		return nullptr;
	}
	inline void printtest(std::ostream & os, const datainfo & info) const {
		test->print(os, info);
	}

#ifdef USEPERSIST
	inline std::array<double, 13> params() const {
		return std::array<double, 13>(
		{stats.n, stats.t, rate,
			stats.sum2, stats.sum[0], stats.sum[1],
			stats.invar(0,0), stats.invar(0,1),
			stats.invar(1,0), stats.invar(1,1),
			mu[0], mu[1], sigma}); 
	}
#else
	inline std::array<double, 8> params() const {
		return std::array<double, 8>(
		{stats.n, stats.t, rate,
			stats.sum2, stats.sum,
			stats.invar,
			mu, sigma}); 
	}
#endif

	inline std::array<double, 2> get_stats() const {
		std::array<double, 2> ret {stats.n, stats.t};
		return ret;
	}

private:
	class ss {
	public:
		double n; // not int, in case need expected value
		double t; // total time
		//  below, y is value of event & x is predictors (1 if not USEPERSIST)
		double sum2; // sum y^2  / n
		pcim::wtT sum;  // sum xy   / n
		pcim::wtvarT invar;  // sum x^2  (==n if not USEPERSIST)
	private:
		friend class boost::serialization::access;
		template<typename Ar>
		void serialize(Ar &ar, const unsigned int ver) {
			if (ver==0) {
				wtT mean;
				ar & BOOST_SERIALIZATION_NVP(n) & BOOST_SERIALIZATION_NVP(t) & BOOST_SERIALIZATION_NVP(mean) & BOOST_SERIALIZATION_NVP(sum2) & BOOST_SERIALIZATION_NVP(invar);
				sum = mean*n;
			} else {
				ar & BOOST_SERIALIZATION_NVP(n) & BOOST_SERIALIZATION_NVP(t) & BOOST_SERIALIZATION_NVP(sum) & BOOST_SERIALIZATION_NVP(sum2) & BOOST_SERIALIZATION_NVP(invar);
			}
		}
	};
	static constexpr size_t nleaffeat=
#ifdef USEPERSIST
		npredfeat+2;
#else
		3;
#endif

#ifdef USEPERSIST
	inline static wtT predfeat(double lastv) {
		wtT ret;
		ret << 1, lastv;
		return ret;
	}
#endif

	pcim(const std::vector<vartrajrange> &data, const pcim::ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params, int m,
		size_t trackbranch = 0);


	void build(const std::vector<vartrajrange> &data, const ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params, size_t trackbranch);
	
	static ss suffstats(const std::vector<vartrajrange> &data);
	static double score(const ss &d, const pcimparams &p);
	void calcleaf(const ss &d, const pcimparams &p);

	double getratevar(const traj &tr, int var, double t, double &until, const pcim *&leaf) const;

	void getstatevar(const traj &tr, int var, double t, double &until,
		size_t & id) const;

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

	size_t leaf_id;// used only for extracting features for coclustering
	shptr<pcim> ftree,ttree;
	shptr<pcimtest> test;
	double rate;
	wtT mu;
	double sigma;
	ss stats;
	int globalm; // number of traj it was built from
#ifdef USEPERSIST
	wtvarT xxinvsqrt;
#endif
public:
	static std::string globalbranch; // branch with most events
	
private:
	void getleaffeature(const std::vector<vartrajrange> &tr,
			std::array<double,nleaffeat> &f) const;

	double llh(const std::vector<vartrajrange> &tr) const;

	void featurenames(std::vector<std::string> &ret, std::string prefix) const;
#ifdef RATEONLY
	std::array<std::string,1> getleaffeaturenames() const;
#else
	std::array<std::string,nleaffeat> getleaffeaturenames() const;
#endif

	void calcxxinvsqrt(const ss &d);

	friend class boost::serialization::access;
/*
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(globalm);
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(mu);
		ar & BOOST_SERIALIZATION_NVP(sigma);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
	}
*/
	template<typename Ar>
	void save(Ar &ar, const unsigned int ver) const {
		ar & BOOST_SERIALIZATION_NVP(globalm);
		ar & BOOST_SERIALIZATION_NVP(globalbranch);
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(mu);
		ar & BOOST_SERIALIZATION_NVP(sigma);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
	}
	template<typename Ar>
	void load(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(globalm);
		ar & BOOST_SERIALIZATION_NVP(globalbranch);
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(mu);
		ar & BOOST_SERIALIZATION_NVP(sigma);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
		calcxxinvsqrt(stats);
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

#ifdef USEPERSIST
namespace boost {
	template<typename Ar>
	inline void serialize(Ar &ar, pcim::wtT &x, const unsigned int v) {
		for(int i=0;i<pcim::npredfeat;i++)
			ar & boost::serialization::make_nvp((std::string("x")+std::to_string(i)).c_str(),x(i));
	}
	template<typename Ar>
	inline void serialize(Ar &ar, pcim::wtvarT &x, const unsigned int v) {
		for(int i=0;i<pcim::npredfeat;i++)
			for(int j=0;j<pcim::npredfeat;j++)
				ar & boost::serialization::make_nvp((std::string("x")+std::to_string(i)+"-"+std::to_string(j)).c_str(),x(i,j));
	}
}
#endif
	

BOOST_CLASS_EXPORT_KEY(pcimtest)
BOOST_CLASS_EXPORT_KEY(lasttest)
BOOST_CLASS_EXPORT_KEY(timetest)
BOOST_CLASS_EXPORT_KEY(counttest)
BOOST_CLASS_EXPORT_KEY(varstattest<counttest>)
BOOST_CLASS_EXPORT_KEY(meantest)
BOOST_CLASS_EXPORT_KEY(varstattest<meantest>)
BOOST_CLASS_EXPORT_KEY(stdvtest)
BOOST_CLASS_EXPORT_KEY(varstattest<stdvtest>)
BOOST_CLASS_EXPORT_KEY(vartest)
BOOST_CLASS_EXPORT_KEY(staticgreqtest)
BOOST_CLASS_EXPORT_KEY(staticeqtest)
BOOST_CLASS_EXPORT_KEY(pcim)
BOOST_SERIALIZATION_SHARED_PTR(pcimtest)
BOOST_SERIALIZATION_SHARED_PTR(pcim)

BOOST_CLASS_VERSION(pcim::ss,1)


#endif
