// next_primes_from_n_fixed.cpp
// Compile: g++ -O3 -std=c++17 next_primes_from_n_fixed.cpp -o next_primes_from_n

#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

using boost::multiprecision::cpp_int;

static const uint64_t SMALL_PRIMES[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
    157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,
    239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,
    331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,
    421,431,433,439,443,449,457,461,463,467,479,487,491,499
};
static const size_t SMALL_PRIME_COUNT = sizeof(SMALL_PRIMES) / sizeof(SMALL_PRIMES[0]);

inline static cpp_int mulmod(const cpp_int& a, const cpp_int& b, const cpp_int& mod) {
  return (a * b) % mod;
}

cpp_int powmod(cpp_int base, cpp_int exp, const cpp_int& mod) {
  cpp_int res = 1 % mod;
  base %= mod;
  while (exp > 0) {
    if ((exp & 1) != 0) res = mulmod(res, base, mod);
    base = mulmod(base, base, mod);
    exp >>= 1;
  }
  return res;
}

static inline void decompose(const cpp_int& n, cpp_int& d, unsigned& s) {
  d = n - 1;
  s = 0;
  while ((d & 1) == 0) { d >>= 1; ++s; }
}

bool miller_rabin(const cpp_int& n, int rounds = 32, std::mt19937_64* rng_ptr = nullptr) {
  if (n < 2) return false;
  for (size_t i = 0; i < SMALL_PRIME_COUNT; ++i) {
    uint64_t p = SMALL_PRIMES[i];
    if (n == p) return true;
    if (n % p == 0) return false;
  }

  cpp_int d; unsigned s;
  decompose(n, d, s);

  std::mt19937_64 local_rng(std::random_device{}());
  std::mt19937_64& rng = rng_ptr ? *rng_ptr : local_rng;
  std::uniform_int_distribution<uint64_t> dist64(2, std::numeric_limits<uint64_t>::max());

  auto random_base = [&](const cpp_int& limit_minus_two)->cpp_int {
    if (limit_minus_two <= 2) return cpp_int(2);
    // calculer le nombre de "limbs" 64-bits nécessaires
    unsigned long msb_index = boost::multiprecision::msb(limit_minus_two);
    unsigned long bits = msb_index + 1;
    size_t limbs = (bits + 63) / 64;
    cpp_int a = 0;
    for (size_t i = 0; i < limbs; ++i) {
      a <<= 64;
      a += cpp_int(dist64(rng));
    }
    // ramener dans [2, limit_minus_two]
    if (a <= 1) a += 2;
    a %= (limit_minus_two - 1);
    a += 2;
    return a;
    };

  for (int t = 0; t < rounds; ++t) {
    cpp_int a = random_base(n - 2);
    cpp_int x = powmod(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool passed = false;
    for (unsigned r = 1; r < s; ++r) {
      x = mulmod(x, x, n);
      if (x == n - 1) { passed = true; break; }
    }
    if (!passed) return false;
  }
  return true;
}

bool is_prime(const cpp_int& n, std::mt19937_64* rng_ptr = nullptr) {
  if (n < 2) return false;
  for (size_t i = 0; i < SMALL_PRIME_COUNT; ++i) {
    uint64_t p = SMALL_PRIMES[i];
    if (n == p) return true;
    if (n % p == 0) return false;
  }
  return miller_rabin(n, 32, rng_ptr);
}

cpp_int next_candidate(cpp_int n) {
  if (n <= 2) return 2;
  if ((n & 1) == 0) ++n;
  while (true) {
    if (n % 3 != 0 && n % 5 != 0) return n;
    n += 2;
  }
}

std::vector<cpp_int> generate_primes(cpp_int start, size_t count) {
  std::vector<cpp_int> primes;
  primes.reserve(count);
  std::mt19937_64 rng(std::random_device{}());
  cpp_int n = next_candidate(start);
  if (n == 2) { primes.push_back(n); n = 3; }
  while (primes.size() < count) {
    if (is_prime(n, &rng)) primes.push_back(n);
    n += 2;
  }
  return primes;
}

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  cpp_int start;
  size_t how_many = 100;

  if (argc >= 2) {
    std::istringstream iss(argv[1]);
    if (!(iss >> start)) {
      std::cerr << "Impossible de lire l'entier de départ.\n";
      return 1;
    }
  }
  else {
    std::istringstream iss("18446744073713598463");
    iss >> start;
  }
  if (argc >= 3) how_many = static_cast<size_t>(std::stoull(argv[2]));

  auto primes = generate_primes(start, how_many);
  for (size_t i = 0; i < primes.size(); ++i) {
    std::cout << primes[i] << '\n';
  }
  return 0;
}
