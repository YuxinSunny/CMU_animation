// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).

template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
                                            const T& position0,
                                            const T& position1,
                                            const T& tangent0,
                                            const T& tangent1,
                                            double normalizedTime,
                                            int derivative )
{
    // TODO IMPLEMENT ME (TASK 1A)
    //*sun*************************

    double t=normalizedTime;
    double h_00,h_11,h_01,h_10,h_00_1,h_00_2,h_01_1,h_01_2,h_10_1,h_10_2,h_11_1,h_11_2,t_2,t_3=0;
    t_2=t*t;
    t_3=t_2*t;
    h_00=2.*t_3-3.*t_2+1;
    h_00_1=6.*t_2-6.*t;
    h_00_2=12.*t-6;
    h_10=t_3-2.*t_2+t;
    h_10_1=3.*t_2-4.*t+1;
    h_10_2=6.*t-4;
    h_01=-2.*t_3+3.*t_2;
    h_01_1=-6.*t_2+6.*t;
    h_01_2=-12.*t+6;
    h_11=t_3-t_2;
    h_11_1=3.*t_2-2.*t;
    h_11_2=6.*t-2;
    if(derivative==0)
        return h_00*position0+h_10*tangent0+h_01*position1+h_11*tangent1;
    else if(derivative==1)
        return h_00_1*position0+h_10_1*tangent0+h_01_1*position1+h_11_1*tangent1;
    else return h_00_2*position0+h_10_2*tangent0+h_01_2*position1+h_11_2*tangent1;
    //*sun*************************
    //return T();//////////sun
}

//// Returns a state interpolated between the values directly before and after the given time.
template <class T>
inline T Spline<T>::evaluate( double time, int derivative )
{
    // TODO IMPLEMENT ME (TASK 1B)

    //*sun********************
    KnotIter it1=knots.begin();
    KnotIter it2=knots.end();
    if(knots.size()<1)
        return T();
    else if(knots.size()==1)
    {
        if(derivative==0)return  it1->second;
        if (derivative == 1||derivative==2) return T();
    }
    else if(time <= it1->first)
    {
        if(derivative==0) return it1->second;
        if (derivative == 1||derivative==2) return T();
    }
    else if(time >= knots.rbegin()->first)//如果直接用end()是指向null的
    {
        if(derivative==0)return knots.rbegin()->second;
        if (derivative == 1||derivative==2) return T();
    }
    else
    {
        double t_0,t_3=0;
        T p_0,p_3;
        KnotIter it_p_1=std::prev(knots.upper_bound(time));
        KnotIter it_p_2=knots.upper_bound(time);
        if(it_p_1==it1)
        {
    
            t_0=2.*it_p_1->first-it_p_2->first;
            p_0=2.*it_p_1->second-it_p_2->second;
        }
        else
        {
            t_0=(std::prev(it_p_1))->first;
            p_0=(std::prev(it_p_1))->second;
        }
        if(it_p_2==(std::prev(it2)))
        {
            t_3=2*it_p_2->first-it_p_1->first;
            p_3=2*it_p_2->second-it_p_1->second;
        }
        else
        {
            t_3=(std::next(it_p_2))->first;
            p_3=(std::next(it_p_2))->second;
        }
    
        T tangent0=(it_p_2->second-p_0)*(it_p_2->first - it_p_1->first)/(it_p_2->first-t_0);
        T tangent1=(p_3-it_p_1->second)*(it_p_2->first - it_p_1->first)/(t_3-it_p_1->first);

        double t=(time-it_p_1->first)/(it_p_2->first-it_p_1->first);
//cuowu
        if (derivative == 1) {
            return cubicSplineUnitInterval(it_p_1->second,it_p_2->second,tangent0,tangent1,t,derivative) / (it_p_2->first - it_p_1->first);
                // return cubicSplineUnitInterval_Centripetal(q0, q1, q2, q3, t, 1) * (1.0 / (t2 - t1));
            }
            else if (derivative == 2) {
                return cubicSplineUnitInterval(it_p_1->second,it_p_2->second,tangent0,tangent1,t,derivative) / pow((it_p_2->first - it_p_1->first), 2);
                            }
            // derivative == 0
        return cubicSplineUnitInterval(it_p_1->second,it_p_2->second,tangent0,tangent1,t,derivative);
    }
    //*sun********************
    //return T();
}




// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance )
{
    // Empty maps have no knots.
    if( knots.size() < 1 )
    {
        return false;
    }
    
    // Look up the first element > or = to time.
    typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
    typename std::map<double, T>::iterator t1_iter;
    t1_iter = t2_iter;
    t1_iter--;
    
    if( t2_iter == knots.end() )
    {
        t2_iter = t1_iter;
    }
    
    // Handle tolerance bounds,
    // because we are working with floating point numbers.
    double t1 = (*t1_iter).first;
    double t2 = (*t2_iter).first;
    
    double d1 = fabs(t1 - time);
    double d2 = fabs(t2 - time);
    
    
    if(d1 < tolerance && d1 < d2)
    {
        knots.erase(t1_iter);
        return true;
    }
    
    if(d2 < tolerance && d2 < d1)
    {
        knots.erase(t2_iter);
        return t2;
    }
    
    return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue( double time, T value )
{
    knots[ time ] = value;
}

template <class T>
inline T Spline<T>::operator()( double time )
{
    return evaluate( time );
}
