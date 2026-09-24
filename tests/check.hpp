#pragma once
#include <cstdio>
#include <cstdlib>
static int g_failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::printf("FEHLER %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); ++g_failures; } } while (0)
#define REPORT() do { if (g_failures) { std::printf("%d Pruefung(en) fehlgeschlagen\n", g_failures); return 1; } std::printf("alle Pruefungen bestanden\n"); return 0; } while (0)
