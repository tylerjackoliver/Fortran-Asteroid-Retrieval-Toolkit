#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

extern "C" {

    void tbpObserver(std::vector<double>& x, const double t){};

    void tbp(std::vector<double>& in, std::vector<double>& out, const double t)
    {

        double d1, d2;
        const double mu = 3.0032080443e-06;

        d1 = sqrt( (in[0]+mu) * (in[0]+mu) + in[1]*in[1] + in[2]*in[2]);
        d1 = 1.0/(d1*d1*d1);    // first distance
        d2 = sqrt( (in[0]+mu-1.0)*(in[0]+mu-1.0) + in[1]*in[1] + in[2]*in[2] );
        d2 = 1.0 / (d2*d2*d2);    // second distance

        out[0] = in[3];
        out[1] = in[4];
        out[2] = in[5]; 
        out[3] = in[0] + 2.0*in[4] - d1*(1-mu)*(in[0]+mu) - d2*mu*(in[0]+mu-1.0);
        out[4] = in[1] - 2.0*in[3] - d1*(1-mu)*in[1]      - d2*mu*in[1];
        out[5] =                   - d1*(1-mu)*in[2]      - d2*mu*in[2];

    }

    void integratorCpp(double input[], const double t)
    {

        typedef std::vector<double> state_type;
        state_type x(6);

        double absTol = 1.e-013;
        double relTol = 1.e-013;

        // Assign values

        x[0] = input[0];
        x[1] = input[1];
        x[2] = input[2];
        x[3] = input[3];
        x[4] = input[4];
        x[5] = input[5];

        // Instantiate stepper

        typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type, double, state_type, double> rk78;
        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled( absTol, relTol, rk78() ), tbp, x, 0., t, -0.0001);

        // Assign back

        input[0] = x[0];
        input[1] = x[1];
        input[2] = x[2];
        input[3] = x[3];
        input[4] = x[4];
        input[5] = x[5];
    
    }

}