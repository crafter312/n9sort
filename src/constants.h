#ifndef CONSTANTS_H
#define CONSTANTS_H

/* This originally used to just be a constants header file, but after modifications
 * by Henry Webb (h.s.webb@wustl.edu) is now more of a general utilities file.
 * Additions include text color codes for cout statements (stolen from TNLIB written
 * by Alex Alafa), hash tables to convert (Z, A) pairs to mass values, and a macro
 * for creating enum classes with automatically generated `ToString` functions.
 */

#include <cmath>
#include <unordered_map>
#include <utility>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

const double m0 = 931.478;
const double c = 30.;
const double vfact = c/sqrt(m0);  // velocity (cm/ns) = vfact * sqrt(2. * E(MeV)/A(amu))
const double pi = acos(-1.);
const double deg_to_rad = pi/180;
const double rad_to_deg = 180/pi;

// Masses excesses from AME2016 compilation
const double mass_n = 8.07132;
const double mass_p = 7.28897;
const double mass_d = 13.13572;
const double mass_t = 14.9498;
const double mass_3He = 14.93121;
const double mass_alpha = 2.42491;
const double mass_5He = 11.231; // AME2020
const double mass_6He = 17.5921;
const double mass_8He = 31.6096;
const double mass_5Li = 11.678886;
const double mass_6Li = 14.0868;
const double mass_7Li = 14.9071;
const double mass_8Li = 20.9458;
const double mass_9Li = 24.9549;
const double mass_6Be = 18.375033;
const double mass_7Be = 15.768999;
const double mass_8Be = 4.9416;
const double mass_9Be = 11.3484;
const double mass_10Be = 12.6074;
const double mass_11Be = 20.1771;
const double mass_8B = 22.9215;
const double mass_9B = 12.416488;
const double mass_10B = 12.0506;
const double mass_11B = 8.6677;
const double mass_8C = 35.064; // AME 2020
const double mass_9C = 28.910972;
const double mass_10C = 15.698672;
const double mass_11C = 10.649396;
const double mass_12C = 0.;
const double mass_13C = 3.12500888;
const double mass_14C = 3.019892;
const double mass_11N = 24.303559;
const double mass_12N = 17.338068;
const double mass_13N = 5.345481;
const double mass_14N = 2.863416;
const double mass_15N = 0.101438;
const double mass_13O = 23.115432;
const double mass_14O = 8.007781;
const double mass_15O = 2.855605;
const double mass_16O = -4.737001;
const double mass_17O = -0.808763;
const double mass_14F = 31.964402;
const double mass_15F = 16.566751;
const double mass_17F = 1.951702;
const double mass_18F = .873113;
const double mass_17Ne = 16.500447;
const double mass_18Ne = 5.317614;

// Total masses
const double Mass_n = m0+mass_n;
const double Mass_p = m0+mass_p;
const double Mass_d = 2.*m0+mass_d;
const double Mass_t = 3.*m0+mass_t;
const double Mass_3He = 3.*m0+mass_3He;
const double Mass_alpha = 4.*m0+mass_alpha;
const double Mass_5He = 5.*m0+mass_5He;
const double Mass_6He = 6.*m0+mass_6He;
const double Mass_8He = 8.*m0+mass_8He;
const double Mass_5Li = 5.*m0+mass_5Li;
const double Mass_6Li = 6.*m0+mass_6Li;
const double Mass_7Li = 7.*m0+mass_7Li;
const double Mass_8Li = 8.*m0+mass_8Li;
const double Mass_9Li = 9.*m0+mass_6Li;
const double Mass_6Be = 6.*m0+mass_6Be;
const double Mass_7Be = 7.*m0+mass_7Be;
const double Mass_8Be = 8.*m0+mass_8Be;
const double Mass_9Be = 9.*m0+mass_9Be;
const double Mass_10Be = 10.*m0+mass_10Be;
const double Mass_11Be = 11.*m0+mass_11Be;
const double Mass_8B = 8.*m0+mass_8B;
const double Mass_9B = 9.*m0+mass_9B;
const double Mass_10B = 10.*m0+mass_10B;
const double Mass_11B = 11.*m0+mass_11B;
const double Mass_8C = 8.*m0+mass_8C;
const double Mass_9C = 9.*m0+mass_9C;
const double Mass_10C = 10.*m0+mass_10C;
const double Mass_11C = 11.*m0+mass_11C;
const double Mass_12C = 12.*m0+mass_12C;
const double Mass_13C = 13.*m0+mass_13C;
const double Mass_14C = 14.*m0+mass_14C;
const double Mass_11N = 11.*m0+mass_11N;
const double Mass_12N = 12.*m0+mass_12N;
const double Mass_13N = 13.*m0+mass_13N;
const double Mass_14N = 14.*m0+mass_14N;
const double Mass_15N = 15.*m0+mass_15N;
const double Mass_13O = 13.*m0+mass_13O;
const double Mass_14O = 14.*m0+mass_14O;
const double Mass_15O = 15.*m0+mass_15O;
const double Mass_16O = 16.*m0+mass_16O;
const double Mass_17O = 17.*m0+mass_17O;
const double Mass_14F = 14.*m0+mass_14F;
const double Mass_15F = 15.*m0+mass_15F;
const double Mass_17F = 17.*m0+mass_17F;
const double Mass_18F = 18.*m0+mass_18F;
const double Mass_17Ne = 17.*m0+mass_17Ne;
const double Mass_18Ne = 18.*m0 + mass_18Ne;

// Pair hasher required for the maps created below
struct pair_hash {
    inline size_t operator()(const std::pair<size_t, size_t> & v) const {
        // A simple way to combine hashes:
        return std::hash<size_t>{}(v.first) ^ (std::hash<size_t>{}(v.second) << 1);
    }
};

// Lookup table for mass excesses (first index Z, second index A)
const std::unordered_map<std::pair<size_t, size_t>, double, pair_hash> mass_lookup = {
	{{0, 1}, mass_n},
	{{1, 1}, mass_p}, {{1, 2}, mass_d}, {{1, 3}, mass_t},
	{{2, 3}, mass_3He}, {{2, 4}, mass_alpha}, {{2, 5}, mass_5He}, {{2, 6}, mass_6He}, {{2, 8}, mass_8He},
	{{3, 5}, mass_5Li}, {{3, 6}, mass_6Li}, {{3, 7}, mass_7Li}, {{3, 8}, mass_8Li}, {{3, 9}, mass_9Li},
	{{4, 6}, mass_6Be}, {{4, 7}, mass_7Be}, {{4, 8}, mass_8Be}, {{4, 9}, mass_9Be}, {{4, 10}, mass_10Be}, {{4, 11}, mass_11Be},
	{{5, 8}, mass_8B}, {{5, 9}, mass_9B}, {{5, 10}, mass_10B}, {{5, 11}, mass_11B},
	{{6, 8}, mass_8C}, {{6, 9}, mass_9C}, {{6, 10}, mass_10C}, {{6, 11}, mass_11C}, {{6, 12}, mass_12C}, {{6, 13}, mass_13C}, {{6, 14}, mass_14C},
	{{7, 11}, mass_11N}, {{7, 12}, mass_12N}, {{7, 13}, mass_13N}, {{7, 14}, mass_14N}, {{7, 15}, mass_15N},
	{{8, 13}, mass_13O}, {{8, 14}, mass_14O}, {{8, 15}, mass_15O}, {{8, 16}, mass_16O}, {{8, 17}, mass_17O},
	{{9, 14}, mass_14F}, {{9, 15}, mass_15F}, {{9, 17}, mass_17F}, {{9, 18}, mass_18F},
	{{10, 17}, mass_17Ne}, {{10, 18}, mass_18Ne}
};

// Lookup table for total masses (first index Z, second index A)
const std::unordered_map<std::pair<size_t, size_t>, double, pair_hash> Mass_lookup = {
	{{0, 1}, Mass_n},
	{{1, 1}, Mass_p}, {{1, 2}, Mass_d}, {{1, 3}, Mass_t},
	{{2, 3}, Mass_3He}, {{2, 4}, Mass_alpha}, {{2, 5}, Mass_5He}, {{2, 6}, Mass_6He}, {{2, 8}, Mass_8He},
	{{3, 5}, Mass_5Li}, {{3, 6}, Mass_6Li}, {{3, 7}, Mass_7Li}, {{3, 8}, Mass_8Li}, {{3, 9}, Mass_9Li},
	{{4, 6}, Mass_6Be}, {{4, 7}, Mass_7Be}, {{4, 8}, Mass_8Be}, {{4, 9}, Mass_9Be}, {{4, 10}, Mass_10Be}, {{4, 11}, Mass_11Be},
	{{5, 8}, Mass_8B}, {{5, 9}, Mass_9B}, {{5, 10}, Mass_10B}, {{5, 11}, Mass_11B},
	{{6, 8}, Mass_8C}, {{6, 9}, Mass_9C}, {{6, 10}, Mass_10C}, {{6, 11}, Mass_11C}, {{6, 12}, Mass_12C}, {{6, 13}, Mass_13C}, {{6, 14}, Mass_14C},
	{{7, 11}, Mass_11N}, {{7, 12}, Mass_12N}, {{7, 13}, Mass_13N}, {{7, 14}, Mass_14N}, {{7, 15}, Mass_15N},
	{{8, 13}, Mass_13O}, {{8, 14}, Mass_14O}, {{8, 15}, Mass_15O}, {{8, 16}, Mass_16O}, {{8, 17}, Mass_17O},
	{{9, 14}, Mass_14F}, {{9, 15}, Mass_15F}, {{9, 17}, Mass_17F}, {{9, 18}, Mass_18F},
	{{10, 17}, Mass_17Ne}, {{10, 18}, Mass_18Ne}
};

// Source - https://stackoverflow.com/a/5094430
// Posted by James McNellis
// Retrieved 2026-04-29, License - CC BY-SA 2.5
/*
#include <boost/preprocessor.hpp>

#define X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE(r, data, elem)    \
    case elem : return BOOST_PP_STRINGIZE(elem);

#define DEFINE_ENUM_WITH_STRING_CONVERSIONS(name, enumerators)                \
    enum name {                                                               \
        BOOST_PP_SEQ_ENUM(enumerators)                                        \
    };                                                                        \
                                                                              \
    static inline const char* ToString(name v)                                       \
    {                                                                         \
        switch (v)                                                            \
        {                                                                     \
            BOOST_PP_SEQ_FOR_EACH(                                            \
                X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE,          \
                name,                                                         \
                enumerators                                                   \
            )                                                                 \
            default: return "[Unknown " BOOST_PP_STRINGIZE(name) "]";         \
        }                                                                     \
    }
*/
// Helper to force values into size_t pairs
template<typename T1, typename T2>
constexpr std::pair<size_t, size_t> sz_pair(T1 z, T2 a) {
    return { static_cast<size_t>(z), static_cast<size_t>(a) };
}

#endif
