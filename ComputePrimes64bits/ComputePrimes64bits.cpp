// next_primes_uint64.cpp
// Trouve `count` nombres premiers >= start (uint64_t).
// Conçu pour MSVC (Visual Studio 2022) et compatible g++/clang.
//
// Usage:
//  - Sous Visual Studio : créer un projet Console, ajouter ce fichier et build/run.
//  - En ligne de commande g++: g++ -O3 -std=c++17 next_primes_uint64.cpp -o next_primes

#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <limits>
#include <string>

// Détection MSVC pour utiliser _umul128
#ifdef _MSC_VER
#include <intrin.h>
#endif

using u64 = uint64_t;
//using u128 = unsigned __int128; // utilisé uniquement sur GCC/Clang

// Multiplie a * b mod m sans overflow pour 64 bits.
// Implémentation portable : utilise _umul128 sur MSVC, __uint128_t sinon.
static inline u64 mul_mod(u64 a, u64 b, u64 mod) {
#ifdef _MSC_VER
  // MSVC : _umul128 retourne la partie basse et la partie haute dans *high
  unsigned long long high;
  unsigned long long low = _umul128(a, b, &high);
  // on a 128-bit = (high << 64) | low. On doit faire (a*b) % mod.
  // Réduction par division 128/64 (utilise builtin 128 via pair high/low).
  // On convertit en builtin 128 pour simplifier la division si possible:
  // MSVC ne supporte __uint128_t ; simulons la division en utilisant
  // l'algorithme de réduction double-and-add — mais pour simplicité et
  // performance, utilisons long double trick si mod fits in 64 bits.
  // Simpler approach: use builtin algorithm: (a*b - q*mod) where q = (unsigned __int128)(a*b)/mod
  // But MSVC lacks unsigned __int128. So we'll implement 128-bit division manually:
  // Use std::uint64_t high, low and perform manual division (Knuth) — but that's long.
  // Simpler and sufficiently fast: use loop doubling (binary method).
  u64 result = 0;
  u64 base = a % mod;
  u64 bb = b;
  while (bb) {
    if (bb & 1) {
      // result = (result + base) % mod, avoid overflow by subtraction
      u64 tmp = result + base;
      if (tmp < result || tmp >= mod) tmp = (tmp % mod);
      result = tmp % mod;
    }

    bb >>= 1;
    if (bb) {
      // base = (base * 2) % mod, safe
      base = base + base;
      if (base >= mod || base < base - base) base %= mod;
    }
  }

  return result % mod;
#else
  // GCC/Clang : on a __uint128_t
  u128 res = (u128)a * (u128)b;
  res %= mod;
  return (u64)res;
#endif
}

// Exponentiation modulaire
static inline u64 pow_mod(u64 a, u64 d, u64 mod) {
  u64 res = 1;
  a %= mod;
  while (d) {
    if (d & 1) res = mul_mod(res, a, mod);
    a = mul_mod(a, a, mod);
    d >>= 1;
  }

  return res;
}

// Miller-Rabin déterministe pour 64-bit (bases spéciales)
static bool is_prime_u64(u64 n) {
  if (n < 2) return false;
  static const u64 small_primes[] = {
      2ull,3ull,5ull,7ull,11ull,13ull,17ull,19ull,23ull,29ull,31ull,37ull
  };
  for (u64 p : small_primes) {
    if (n == p) return true;
    if (n % p == 0) return false;
  }

  // écrire n-1 = d * 2^s
  u64 d = n - 1;
  int s = 0;
  while ((d & 1) == 0) { d >>= 1; ++s; }

  // bases déterministes pour n < 2^64
  const u64 bases[] = { 2ull, 325ull, 9375ull, 28178ull, 450775ull, 9780504ull, 1795265022ull };

  for (u64 a : bases) {
    if (a % n == 0) continue;
    u64 x = pow_mod(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool composite = true;
    for (int r = 1; r < s; ++r) {
      x = mul_mod(x, x, n);
      if (x == n - 1) { composite = false; break; }
    }
    if (composite) return false;
  }

  return true;
}

// Renvoie le prochain candidat impair >= n
static u64 next_candidate_u64(u64 n) {
  if (n <= 2) return 2;
  if ((n & 1) == 0) ++n;
  // skip multiples of 3 quickly
  while (n % 3 == 0) n += 2;
  return n;
}

static std::vector<u64> generate_primes_u64(u64 start, size_t count) {
  std::vector<u64> primes;
  primes.reserve(count);
  u64 n = next_candidate_u64(start);
  if (n == 2) { primes.push_back(2); n = 3; }
  while (primes.size() < count) {
    if (is_prime_u64(n)) primes.push_back(n);
    if (n >= std::numeric_limits<u64>::max() - 2) break; // sécurité
    n += 2;
  }

  return primes;
}

int main(int argc, char** argv) {
  u64 start = 18446744073709551615ULL; // exemple fourni
  size_t count = 100;
  if (argc >= 2) {
    // lire en decimal (potentiellement grand)
    try {
      std::string s = argv[1];
      // stoi pour unsigned long long
      start = std::stoull(s);
    }
    catch (...) {
      std::cerr << "Argument invalide pour start\n";
      return 1;
    }
  }
  if (argc >= 3) {
    count = static_cast<size_t>(std::stoull(argv[2]));
  }

  auto primes = generate_primes_u64(start, count);
  for (u64 p : primes) std::cout << p << '\n';
  return 0;
}
