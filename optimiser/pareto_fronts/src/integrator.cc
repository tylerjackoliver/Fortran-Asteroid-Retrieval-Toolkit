#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

extern "C" {

    /* @brief Spoof observer function for boost; does nothing */
    void tbpObserver(std::vector<double>& x, const double t){};

    /* @brief Force-function for the CR3BP
     * @param[in] in std::vector<double> of current state vector
     * @param[inout] out std::vector<double> of derivative of current state vector
     * @param[in] t Current time
     */
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

    /* @brief Integrate a given state backwards in the CR3BP
     * @param[inout] input Double-valued array of values to integrate
     * @param[in] t Time to integrate backwards to
     */
    void integratorCpp(double input[], const double t)
    {
        // Define state for the integration
        typedef std::vector<double> state_type;
        state_type x(6);

        // Define tolerances
        double absTol = 1.e-013;
        double relTol = 1.e-013;

        // Assign values
        for (size_t idx = 0; idx < x.size(); ++idx) x[idx] = input[idx];

        // Instantiate stepper type
        typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type, double, state_type, double> rk78;
        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled( absTol, relTol, rk78() ), tbp, x, 0., t, -0.0001);

        // Assign back
        for (size_t idx = 0; idx < x.size(); ++idx) input[idx] = x[idx];
    }
}