#include "parType.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

parType::parType(int Z0, int A0, detType det) : Z(Z0), A(A0), detector(det) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void parType::zeroMask() {
  for (size_t i = 0; i < 6; i++) mask[i] = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void parType::setMask() {
  for (size_t i = 0; i < 6; i++) mask[i] = true;
}



