
// Code generated by stanc f556d0df
#include <stan/model/model_header.hpp>
namespace ppca_model_namespace {

template <typename T, typename S>
std::vector<T> resize_to_match__(std::vector<T>& dst, const std::vector<S>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, -1>
resize_to_match__(Eigen::Matrix<T, -1, -1>& dst, const Eigen::Matrix<T, -1, -1>& src) {
  dst.resize(src.rows(), src.cols());
  return dst;
}

template <typename T>
Eigen::Matrix<T, 1, -1>
resize_to_match__(Eigen::Matrix<T, 1, -1>& dst, const Eigen::Matrix<T, 1, -1>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, 1>
resize_to_match__(Eigen::Matrix<T, -1, 1>& dst, const Eigen::Matrix<T, -1, 1>& src) {
  dst.resize(src.size());
  return dst;
}
std::vector<double> to_doubles__(std::initializer_list<double> x) {
  return x;
}

std::vector<stan::math::var> to_vars__(std::initializer_list<stan::math::var> x) {
  return x;
}

inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}

inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}


using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math; 

static int current_statement__ = 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 9, column 2 to column 17)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 11, column 2 to column 17)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 13, column 2 to column 22)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 15, column 2 to column 15)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 17, column 2 to column 28)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 21, column 2 to column 30)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 23, column 4 to column 36)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 22, column 16 to line 24, column 3)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 22, column 2 to line 24, column 3)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 25, column 2 to column 26)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 26, column 2 to column 26)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 30, column 4 to column 45)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 29, column 16 to line 31, column 3)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 29, column 2 to line 31, column 3)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 2, column 2 to column 19)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 3, column 2 to column 19)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 4, column 2 to column 19)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/stan_compute/ppca.stan', line 5, column 2 to column 18)"};



class ppca_model : public model_base_crtp<ppca_model> {

 private:
  int pos__;
  int N;
  int D;
  int M;
  std::vector<Eigen::Matrix<double, -1, 1>> x;
 
 public:
  ~ppca_model() { }
  
  std::string model_name() const { return "ppca_model"; }
  
  ppca_model(stan::io::var_context& context__,
             unsigned int random_seed__ = 0,
             std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    typedef double local_scalar_t__;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "ppca_model_namespace::ppca_model";
    (void) function__;  // suppress unused var warning
    
    try {
      
      pos__ = 1;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      
      current_statement__ = 15;
      N = context__.vals_i("N")[(1 - 1)];
      context__.validate_dims("data initialization","D","int",
          context__.to_vec());
      
      current_statement__ = 16;
      D = context__.vals_i("D")[(1 - 1)];
      context__.validate_dims("data initialization","M","int",
          context__.to_vec());
      
      current_statement__ = 17;
      M = context__.vals_i("M")[(1 - 1)];
      current_statement__ = 18;
      validate_non_negative_index("x", "N", N);
      current_statement__ = 18;
      validate_non_negative_index("x", "D", D);
      context__.validate_dims("data initialization","x","double",
          context__.to_vec(N, D));
      x = std::vector<Eigen::Matrix<double, -1, 1>>(N, Eigen::Matrix<double, -1, 1>(D));
      
      {
        std::vector<local_scalar_t__> x_flat__;
        current_statement__ = 18;
        assign(x_flat__, nil_index_list(), context__.vals_r("x"),
          "assigning variable x_flat__");
        current_statement__ = 18;
        pos__ = 1;
        current_statement__ = 18;
        for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
          current_statement__ = 18;
          for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 18;
            assign(x,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              x_flat__[(pos__ - 1)], "assigning variable x");
            current_statement__ = 18;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 15;
      current_statement__ = 15;
      check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 16;
      current_statement__ = 16;
      check_greater_or_equal(function__, "D", D, 0);
      current_statement__ = 17;
      current_statement__ = 17;
      check_greater_or_equal(function__, "M", M, 0);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      current_statement__ = 1;
      validate_non_negative_index("z", "M", M);
      current_statement__ = 1;
      validate_non_negative_index("z", "N", N);
      num_params_r__ += M * N;
      current_statement__ = 2;
      validate_non_negative_index("w", "D", D);
      current_statement__ = 2;
      validate_non_negative_index("w", "M", M);
      num_params_r__ += D * M;
      num_params_r__ += 1;
      current_statement__ = 4;
      validate_non_negative_index("mu", "D", D);
      num_params_r__ += D;
      current_statement__ = 5;
      validate_non_negative_index("alpha", "M", M);
      num_params_r__ += M;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename T__>
  T__ log_prob(std::vector<T__>& params_r__, std::vector<int>& params_i__,
               std::ostream* pstream__ = 0) const {
    typedef T__ local_scalar_t__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "ppca_model_namespace::log_prob";
(void) function__;  // suppress unused var warning

    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    
    try {
      current_statement__ = 1;
      validate_non_negative_index("z", "M", M);
      current_statement__ = 1;
      validate_non_negative_index("z", "N", N);
      Eigen::Matrix<local_scalar_t__, -1, -1> z;
      z = Eigen::Matrix<local_scalar_t__, -1, -1>(M, N);
      
      current_statement__ = 1;
      z = in__.matrix(M, N);
      current_statement__ = 2;
      validate_non_negative_index("w", "D", D);
      current_statement__ = 2;
      validate_non_negative_index("w", "M", M);
      Eigen::Matrix<local_scalar_t__, -1, -1> w;
      w = Eigen::Matrix<local_scalar_t__, -1, -1>(D, M);
      
      current_statement__ = 2;
      w = in__.matrix(D, M);
      local_scalar_t__ sigma;
      
      current_statement__ = 3;
      sigma = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        sigma = stan::math::lb_constrain(sigma, 0, lp__);
      } else {
        current_statement__ = 3;
        sigma = stan::math::lb_constrain(sigma, 0);
      }
      current_statement__ = 4;
      validate_non_negative_index("mu", "D", D);
      Eigen::Matrix<local_scalar_t__, -1, 1> mu;
      mu = Eigen::Matrix<local_scalar_t__, -1, 1>(D);
      
      current_statement__ = 4;
      mu = in__.vector(D);
      current_statement__ = 5;
      validate_non_negative_index("alpha", "M", M);
      Eigen::Matrix<local_scalar_t__, -1, 1> alpha;
      alpha = Eigen::Matrix<local_scalar_t__, -1, 1>(M);
      
      current_statement__ = 5;
      alpha = in__.vector(M);
      current_statement__ = 5;
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        current_statement__ = 5;
        if (jacobian__) {
          current_statement__ = 5;
          assign(alpha, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(alpha[(sym1__ - 1)], 0, lp__),
            "assigning variable alpha");
        } else {
          current_statement__ = 5;
          assign(alpha, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(alpha[(sym1__ - 1)], 0),
            "assigning variable alpha");
        }}
      {
        current_statement__ = 6;
        lp_accum__.add(normal_log<propto__>(to_vector(z), 0, 1));
        current_statement__ = 9;
        for (size_t d = 1; d <= D; ++d) {
          current_statement__ = 7;
          lp_accum__.add(
            normal_log<propto__>(
              rvalue(w, cons_list(index_uni(d), nil_index_list()), "w"), 0,
              multiply(sigma, alpha)));}
        current_statement__ = 10;
        lp_accum__.add(lognormal_log<propto__>(sigma, 0, 1));
        current_statement__ = 11;
        lp_accum__.add(inv_gamma_log<propto__>(alpha, 1, 1));
        current_statement__ = 14;
        for (size_t n = 1; n <= N; ++n) {
          current_statement__ = 12;
          lp_accum__.add(
            normal_log<propto__>(x[(n - 1)], add(multiply(w, col(z, n)), mu),
              sigma));}
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob() 
    
  template <typename RNG>
  void write_array(RNG& base_rng__, std::vector<double>& params_r__,
                   std::vector<int>& params_i__, std::vector<double>& vars__,
                   bool emit_transformed_parameters__ = true,
                   bool emit_generated_quantities__ = true,
                   std::ostream* pstream__ = 0) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "ppca_model_namespace::write_array";
(void) function__;  // suppress unused var warning

    (void) function__;  // suppress unused var warning

    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    
    try {
      current_statement__ = 1;
      validate_non_negative_index("z", "M", M);
      current_statement__ = 1;
      validate_non_negative_index("z", "N", N);
      Eigen::Matrix<double, -1, -1> z;
      z = Eigen::Matrix<double, -1, -1>(M, N);
      
      current_statement__ = 1;
      z = in__.matrix(M, N);
      current_statement__ = 2;
      validate_non_negative_index("w", "D", D);
      current_statement__ = 2;
      validate_non_negative_index("w", "M", M);
      Eigen::Matrix<double, -1, -1> w;
      w = Eigen::Matrix<double, -1, -1>(D, M);
      
      current_statement__ = 2;
      w = in__.matrix(D, M);
      double sigma;
      
      current_statement__ = 3;
      sigma = in__.scalar();
      current_statement__ = 3;
      sigma = stan::math::lb_constrain(sigma, 0);
      current_statement__ = 4;
      validate_non_negative_index("mu", "D", D);
      Eigen::Matrix<double, -1, 1> mu;
      mu = Eigen::Matrix<double, -1, 1>(D);
      
      current_statement__ = 4;
      mu = in__.vector(D);
      current_statement__ = 5;
      validate_non_negative_index("alpha", "M", M);
      Eigen::Matrix<double, -1, 1> alpha;
      alpha = Eigen::Matrix<double, -1, 1>(M);
      
      current_statement__ = 5;
      alpha = in__.vector(M);
      current_statement__ = 5;
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        current_statement__ = 5;
        assign(alpha, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(alpha[(sym1__ - 1)], 0),
          "assigning variable alpha");}
      for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= M; ++sym2__) {
          vars__.push_back(
            rvalue(z,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())), "z"));}}
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= D; ++sym2__) {
          vars__.push_back(
            rvalue(w,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())), "w"));}}
      vars__.push_back(sigma);
      for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
        vars__.push_back(mu[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        vars__.push_back(alpha[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array() 
    
  void transform_inits(const stan::io::var_context& context__,
                       std::vector<int>& params_i__,
                       std::vector<double>& vars__, std::ostream* pstream__) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      
      pos__ = 1;
      current_statement__ = 1;
      validate_non_negative_index("z", "M", M);
      current_statement__ = 1;
      validate_non_negative_index("z", "N", N);
      Eigen::Matrix<double, -1, -1> z;
      z = Eigen::Matrix<double, -1, -1>(M, N);
      
      {
        std::vector<local_scalar_t__> z_flat__;
        current_statement__ = 1;
        assign(z_flat__, nil_index_list(), context__.vals_r("z"),
          "assigning variable z_flat__");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 1;
          for (size_t sym2__ = 1; sym2__ <= M; ++sym2__) {
            current_statement__ = 1;
            assign(z,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              z_flat__[(pos__ - 1)], "assigning variable z");
            current_statement__ = 1;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 2;
      validate_non_negative_index("w", "D", D);
      current_statement__ = 2;
      validate_non_negative_index("w", "M", M);
      Eigen::Matrix<double, -1, -1> w;
      w = Eigen::Matrix<double, -1, -1>(D, M);
      
      {
        std::vector<local_scalar_t__> w_flat__;
        current_statement__ = 2;
        assign(w_flat__, nil_index_list(), context__.vals_r("w"),
          "assigning variable w_flat__");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
          current_statement__ = 2;
          for (size_t sym2__ = 1; sym2__ <= D; ++sym2__) {
            current_statement__ = 2;
            assign(w,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              w_flat__[(pos__ - 1)], "assigning variable w");
            current_statement__ = 2;
            pos__ = (pos__ + 1);}}
      }
      double sigma;
      
      current_statement__ = 3;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      current_statement__ = 3;
      sigma = stan::math::lb_free(sigma, 0);
      current_statement__ = 4;
      validate_non_negative_index("mu", "D", D);
      Eigen::Matrix<double, -1, 1> mu;
      mu = Eigen::Matrix<double, -1, 1>(D);
      
      {
        std::vector<local_scalar_t__> mu_flat__;
        current_statement__ = 4;
        assign(mu_flat__, nil_index_list(), context__.vals_r("mu"),
          "assigning variable mu_flat__");
        current_statement__ = 4;
        pos__ = 1;
        current_statement__ = 4;
        for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
          current_statement__ = 4;
          assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
            mu_flat__[(pos__ - 1)], "assigning variable mu");
          current_statement__ = 4;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 5;
      validate_non_negative_index("alpha", "M", M);
      Eigen::Matrix<double, -1, 1> alpha;
      alpha = Eigen::Matrix<double, -1, 1>(M);
      
      {
        std::vector<local_scalar_t__> alpha_flat__;
        current_statement__ = 5;
        assign(alpha_flat__, nil_index_list(), context__.vals_r("alpha"),
          "assigning variable alpha_flat__");
        current_statement__ = 5;
        pos__ = 1;
        current_statement__ = 5;
        for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
          current_statement__ = 5;
          assign(alpha, cons_list(index_uni(sym1__), nil_index_list()),
            alpha_flat__[(pos__ - 1)], "assigning variable alpha");
          current_statement__ = 5;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 5;
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        current_statement__ = 5;
        assign(alpha, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(alpha[(sym1__ - 1)], 0),
          "assigning variable alpha");}
      for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= M; ++sym2__) {
          vars__.push_back(
            rvalue(z,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())), "z"));}}
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= D; ++sym2__) {
          vars__.push_back(
            rvalue(w,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())), "w"));}}
      vars__.push_back(sigma);
      for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
        vars__.push_back(mu[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
        vars__.push_back(alpha[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits() 
    
  void get_param_names(std::vector<std::string>& names__) const {
    
    names__.resize(0);
    names__.push_back("z");
    names__.push_back("w");
    names__.push_back("sigma");
    names__.push_back("mu");
    names__.push_back("alpha");
    } // get_param_names() 
    
  void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.resize(0);
    std::vector<size_t> dims__;
    dims__.push_back(M);
    
    dims__.push_back(N);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(D);
    
    dims__.push_back(M);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(D);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(M);
    dimss__.push_back(dims__);
    dims__.resize(0);
    
    } // get_dims() 
    
  void constrained_param_names(std::vector<std::string>& param_names__,
                               bool emit_transformed_parameters__ = true,
                               bool emit_generated_quantities__ = true) const {
    
    for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
      {
        for (size_t sym2__ = 1; sym2__ <= M; ++sym2__) {
          {
            param_names__.push_back(std::string() + "z" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
          }}
      }}
    for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
      {
        for (size_t sym2__ = 1; sym2__ <= D; ++sym2__) {
          {
            param_names__.push_back(std::string() + "w" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
          }}
      }}
    param_names__.push_back(std::string() + "sigma");
    for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
      {
        param_names__.push_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
      {
        param_names__.push_back(std::string() + "alpha" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  void unconstrained_param_names(std::vector<std::string>& param_names__,
                                 bool emit_transformed_parameters__ = true,
                                 bool emit_generated_quantities__ = true) const {
    
    for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
      {
        for (size_t sym2__ = 1; sym2__ <= M; ++sym2__) {
          {
            param_names__.push_back(std::string() + "z" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
          }}
      }}
    for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
      {
        for (size_t sym2__ = 1; sym2__ <= D; ++sym2__) {
          {
            param_names__.push_back(std::string() + "w" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
          }}
      }}
    param_names__.push_back(std::string() + "sigma");
    for (size_t sym1__ = 1; sym1__ <= D; ++sym1__) {
      {
        param_names__.push_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    for (size_t sym1__ = 1; sym1__ <= M; ++sym1__) {
      {
        param_names__.push_back(std::string() + "alpha" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"z\",\"type\":{\"name\":\"matrix\",\"rows\":" << M << ",\"cols\":" << N << "},\"block\":\"parameters\"},{\"name\":\"w\",\"type\":{\"name\":\"matrix\",\"rows\":" << D << ",\"cols\":" << M << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"parameters\"},{\"name\":\"alpha\",\"type\":{\"name\":\"vector\",\"length\":" << M << "},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"z\",\"type\":{\"name\":\"matrix\",\"rows\":" << M << ",\"cols\":" << N << "},\"block\":\"parameters\"},{\"name\":\"w\",\"type\":{\"name\":\"matrix\",\"rows\":" << D << ",\"cols\":" << M << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"parameters\"},{\"name\":\"alpha\",\"type\":{\"name\":\"vector\",\"length\":" << M << "},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool emit_transformed_parameters__ = true,
                     bool emit_generated_quantities__ = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng__, params_r_vec, params_i_vec, vars_vec,
          emit_transformed_parameters__, emit_generated_quantities__, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    template <bool propto__, bool jacobian__, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto__,jacobian__,T_>(vec_params_r, vec_params_i, pstream);
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }

};
}
typedef ppca_model_namespace::ppca_model stan_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif

