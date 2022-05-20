#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "functions.h"

using namespace std;

// TODO Write this function
int energy(const Pixel *const*image, int col, int row, int width, int height)
{
  /*
  The energy of pixel (x, y) is Δx2(x, y) + Δy2(x, y), 
  where the square of the x-gradient Δx2(x, y) = Rx(x, y)2 + Gx(x, y)2 + Bx(x, y)2, 
  and where the central differences Rx(x, y), Gx(x, y), and Bx(x, y) are the absolute value in 
  differences of red, green, and blue components between pixel (x   + 1, y) and pixel (x− 1, y). 
  */
  int energy;

  int red_x; // = abs((image[col + 1][row].r) - (image[col - 1][row].r));
  int blue_x; // = abs((image[col + 1][row].b) - (image[col - 1][row].b));
  int green_x; // = abs((image[col + 1][row].g) - (image[col - 1][row].g));
  
  int red_y; // = abs((image[col][row + 1].r) - (image[col][row - 1].r));
  int blue_y; // = abs((image[col][row + 1].b) - (image[col][row - 1].b));
  int green_y; // = abs((image[col][row + 1].g) - (image[col][row - 1].g));

  if (width > 1) {
    if (col == width - 1) {
      red_x = abs((image[0][row].r) - (image[col - 1][row].r));
      green_x = abs((image[0][row].g) - (image[col - 1][row].g));
      blue_x = abs((image[0][row].b) - (image[col - 1][row].b));
    }
    else if (col == 0) {
      red_x = abs((image[col + 1][row].r) - (image[width - 1][row].r));
      green_x = abs((image[col + 1][row].g) - (image[width - 1][row].g));
      blue_x = abs((image[col + 1][row].b) - (image[width - 1][row].b));
    }
    else {
      red_x = abs((image[col + 1][row].r) - (image[col - 1][row].r));
      green_x = abs((image[col + 1][row].g) - (image[col - 1][row].g));
      blue_x = abs((image[col + 1][row].b) - (image[col - 1][row].b));
    }
  }
  else {
    red_x = 0;
    green_x = 0;
    blue_x = 0;
  }

  if (height > 1) {
    if (row == height - 1) {
      red_y = abs((image[col][0].r) - (image[col][row - 1].r));
      green_y = abs((image[col][0].g) - (image[col][row - 1].g));
      blue_y = abs((image[col][0].b) - (image[col][row - 1].b));
    }
    else if (row == 0) {
      red_y = abs((image[col][row + 1].r) - (image[col][height - 1].r));
      green_y = abs((image[col][row + 1].g) - (image[col][height - 1].g));
      blue_y = abs((image[col][row + 1].b) - (image[col][height - 1].b));
    }
    else {
      red_y = abs((image[col][row + 1].r) - (image[col][row - 1].r));
      green_y = abs((image[col][row + 1].g) - (image[col][row - 1].g));
      blue_y = abs((image[col][row + 1].b) - (image[col][row - 1].b));
    }
  }
  else {
    red_y = 0;
    green_y = 0;
    blue_y = 0;
  }
  int x = pow(red_x, 2) + pow(green_x, 2) + pow(blue_x, 2);
  int y = pow(red_y, 2) + pow(green_y, 2) + pow(blue_y, 2);

  energy = x + y;
  return energy;
}

// TODO Write this function
int getVerticalSeam(const Pixel *const*image, int start_col, int width, int height, int* seam)
{
  
  int row = 0;
  int collumn = start_col;
  int current  = energy(image, collumn, row, width, height);
  int total_energy = current;
  int left;
  int right;
  seam[row] = collumn;
  //std::cout << "This loop should run " << height << " times." << endl; // debug
  //std::cout << row << ": " << "Collumn: " << collumn << " Energy: " << total_energy << endl; //debug
    for (row = 1; row < height; row++) {
      current = (energy(image, collumn, row, width, height)); 
      //cout << row << ": collumn: " << collumn << " move forward at "; //debug
      if (collumn == width - 1) {
        right = energy(image, (collumn - 1), row, width, height);
        left = current;
        //cout << current << " or " << right; // debug
      }
      else if (collumn == 0) {
        left = energy(image, (collumn + 1), row, width, height);
        right = left;
        //cout << left << " or " << current; //debug
      }
      else {
        right = energy(image, (collumn - 1), row, width, height);
        left = energy(image, (collumn + 1), row, width, height);
        //cout << left << " or " << current << " or " << right; // debug
      }
      // move places and add energy
      if (current <= left && current <= right) {
        seam[row] = collumn;
        total_energy += current;
        //cout << " chose " << current << endl; //debug
      }
      else if (left < current && left <= right) {
        current = left;
        collumn += 1;
        seam[row] = collumn;
        total_energy += current;
        //cout << " chose " << left << endl; //debug
      }
      else {
        current = right;
        collumn -= 1;
        seam[row] = collumn;
        total_energy += current;
        //cout << " chose " << right << endl; //debug
      }
      //cout << "new collumn: " << collumn << " new = " << total_energy << endl;
    }

  

  //std::cout << "TOTAL: " << total_energy << endl; //debug

  return total_energy;
  
 /*
 int row = 0;
 int total_energy = energy(image, start_col, row, width, height);
 int collumn = start_col;
 int ahead = 0;

 while (row < height - 1) {
  //int left = energy(image, collumn + 1, row, width, height);
  //int right = energy(image, collumn - 1, row, width, height);
  ahead = energy(image, collumn, row + 1, width, height);
  if (collumn < width - 1 && collumn > 0) {
    int left = energy(image, collumn + 1, row + 1, width, height);
    int right = energy(image, collumn - 1, row + 1, width, height);
    if (left < ahead && left <= right) {
      left = ahead;
      collumn += 1;
    }
    else if (right < left && right < ahead) {
      right = ahead;
      collumn -= 1;
    }
  }
  else if (collumn == width - 1) {
    int right = energy(image, collumn - 1, row + 1, width, height);
    if (right < ahead) {
      right = ahead;
      collumn -= 1;
    }
  }
  else {
    int left = energy(image, collumn + 1, row + 1, width, height);
    if (left < ahead) {
      left = ahead;
      collumn += 1;
    }
  }
  total_energy += ahead;
  row++;
 }
 return total_energy;
 */
}

// TODO Write this function
void removeVerticalSeam(Pixel **image, int& width, int height, int *verticalSeam)
{
  int row = 0;
  for (row = 0; row < height; row++) {
    int seam;
    seam = verticalSeam[row];
    Pixel value;
    value = image[seam][row];
    for (int col = seam; col < width; col++) {
      if (col != width - 1) {
        image[col][row] = image[col + 1][row];
      }
    }
    image[width - 1][row] = value;
  }
  width --;
}

// TODO Write this function for extra credit
int getHorizontalSeam(const Pixel *const*image, int start_row, int width, int height, int* seam)
{
  int col = 0;
  int row = start_row;
  int current  = energy(image, col, row, width, height);
  int total_energy = current;
  int up;
  int down;
  seam[col] = row;
  //std::cout << "This loop should run " << height << " times." << endl; // debug
  //std::cout << row << ": " << "Collumn: " << collumn << " Energy: " << total_energy << endl; //debug
  for (col = 1; col < width; col++) {
    current = (energy(image, col, row, width, height)); 
    //cout << row << ": collumn: " << collumn << " move forward at "; //debug
    if (row == height - 1) {
      down = energy(image, col, (row - 1), width, height);
      up = current;
      //cout << current << " or " << right; // debug
    }
    else if (row == 0) {
      up = energy(image, col, (row + 1), width, height);
      down = up;
      //cout << left << " or " << current; //debug
    }
    else {
      down = energy(image, col, (row - 1), width, height);
      up = energy(image, col, (row + 1), width, height);
      //cout << left << " or " << current << " or " << right; // debug
    }
    // move places and add energy
    if (current <= up && current <= down) {
      seam[col] = row;
      total_energy += current;
      //cout << " chose " << current << endl; //debug
    }
    else if (up < current && up <= down) {
      current = up;
      row += 1;
      seam[col] = row;
      total_energy += current;
      //cout << " chose " << left << endl; //debug
    }
    else {
      current = down;
      row -= 1;
      seam[col] = row;
      total_energy += current;
      //cout << " chose " << right << endl; //debug
    }
    //cout << "new collumn: " << collumn << " new = " << total_energy << endl;
  }
  //std::cout << "TOTAL: " << total_energy << endl; //debug

  return total_energy;
}


// TODO Write this function for extra credit
void removeHorizontalSeam(Pixel **image, int width, int& height, int *horizontalSeam)
{
  int col = 0;
  for (col = 0; col < width; col++) {
    int seam;
    seam = horizontalSeam[col];
    Pixel value;
    value = image[col][seam];
    for (int row = seam; row < height; row++) {
      if (row != height - 1) {
        image[col][row] = image[col][row + 1];
      }
    }
    image[col][height - 1] = value;
  }
  height --;
}

int *findMinVerticalSeam(const Pixel *const*image, int width, int height)
{
  // initialize minSeam and minDistance to seam starting at first col (index 0)
  int *minSeam = new int[height]{0};
  int minDist = getVerticalSeam(image, 0, width, height, minSeam);

  int *candidateSeam = new int[height]{0};
  int candidateDistance = -1; // invalid distance

  // start at second col (index 1) since we initialized with first col (index 0)
  for (int col = 1; col < width; ++col)
  {
    candidateDistance = getVerticalSeam(image, col, width, height, candidateSeam);

    if (candidateDistance < minDist)
    { // new min
      //  swap min & candidate
      minDist = candidateDistance;
      int* temp = candidateSeam;
      candidateSeam = minSeam;
      minSeam = temp;
    }
  }

  // clean up 
  delete [] candidateSeam;

  return minSeam;
}

int *findMinHorizontalSeam(const Pixel *const*image, int width, int height)
{
  // initialize minSeam and minDistance to seam starting at first row (index 0)
  int *minSeam = new int[width]{0};
  int minDistance = getHorizontalSeam(image, 0, width, height, minSeam);

  int *candidateSeam = new int[width]{0};
  int candidateDistance = -1; // invalid distance

  // start at second row (index 1) since we initialized with first row (index 0)
  for (int row = 1; row < height; ++row)
  {
    candidateDistance = getHorizontalSeam(image, row, width, height, candidateSeam);

    if (candidateDistance < minDistance)
    { // new minimum
      //  swap min and candidate seams
      minDistance = candidateDistance;
      int* temp = candidateSeam;
      candidateSeam = minSeam;
      minSeam = temp;
    }
  }

    // clean up 
  delete [] candidateSeam;

  return minSeam;
}

Pixel **createImage(int width, int height)
{
  cout << "Start createImage... " << endl;

  // Create a one dimensional array on the heap of pointers to Pixels
  //    that has width elements (i.e. the number of columns)
  Pixel **image = new Pixel *[width] {}; // initializes to nullptr

  for (int col = 0; col < width; ++col)
  { // loop through each column
    // assign that column to a one dimensional array on the heap of Pixels
    //  that has height elements (i.e. the number of rows)
    try
    {
      image[col] = new Pixel[height];
    }
    catch (std::bad_alloc &e)
    {
      // clean up already allocated arrays
      for (int i = 0; i < col; ++i)
      {
        delete[] image[i];
      }
      delete[] image;
      // rethrow
      throw e;
    }
  }

  // initialize cells
  // cout << "Initializing cells..." << endl;
  for (int row = 0; row < height; ++row)
  {
    for (int col = 0; col < width; ++col)
    {
      // cout << "(" << col << ", " << row << ")" << endl;
      image[col][row] = {0, 0, 0};
    }
  }
  cout << "End createImage... " << endl;
  return image;
}

void deleteImage(Pixel **image, int width)
{
  cout << "Start deleteImage..." << endl;
  // avoid memory leak by deleting the array
  for (int i = 0; i < width; ++i)
  {
    delete[] image[i];
  }
  delete[] image;
  image = nullptr;
  cout << "End deleteImage..." << endl;
}

bool isValidColor(int colorVal)
{
  if (colorVal < 0 || colorVal > 255)
  {
    return false;
  }
  return true;
}

Pixel ** loadImage(string filename, int &width, int &height)
{
  cout << "Start loadImage..." << endl;
  // remove
  ifstream ifs(filename);
  if (!ifs.is_open())
  {
    throw std::invalid_argument("Failed to open input file (" + filename + ")");
  }

  string type;
  ifs >> type; // should be P3
  if (toupper(type.at(0)) != 'P' || type.at(1) != '3')
  {
    throw std::domain_error("Not PPM type P3 (" + type + ")");
  }
  ifs >> width;
  // cout << "w and h: " << w << " " << h << endl;
  if (ifs.fail())
  {
    throw std::domain_error("Read non-integer value for width");
  }
  if (width <= 0)
  {
    ostringstream oss;
    oss << "Width in file must be greater than 0 (" << width << ")";
    throw std::domain_error(oss.str());
  }

  ifs >> height;
  if (ifs.fail())
  {
    cout << "Read non-integer value for height" << endl;
  }
  if (height <= 0)
  {
    ostringstream oss;
    oss << "Height in file must be greater than 0 (" << height << ")";
    throw std::domain_error(oss.str());
  }

  int colorMax = 0;
  ifs >> colorMax;
  if (ifs.fail())
  {
    throw std::domain_error("Read non-integer value for max color value");
  }
  if (colorMax != 255)
  {
    ostringstream oss;
    oss << "Max color value must be 255 (" << colorMax << ")";
    throw std::domain_error(oss.str());
  }

  // load image throws exceptions but we will let them pass through
  Pixel **image = createImage(width, height);

  for (int row = 0; row < height; ++row)
  {
    for (int col = 0; col < width; ++col)
    {
      // cout << "Pixel(" << col << ", " << row << ")" << endl;
      ifs >> image[col][row].r;
      if (ifs.fail() && !ifs.eof())
      { // bad input that is not caused by being at the end of the file
        throw std::domain_error("Read non-integer value for red");
      }
      if (!isValidColor(image[col][row].r))
      {
        ostringstream oss;
        oss << "Invalid color value for red (" << image[col][row].r << ")";
        throw std::domain_error(oss.str());
      }

      ifs >> image[col][row].g;
      if (ifs.fail() && !ifs.eof())
      { // bad input that is not caused by being at the end of the file
        throw std::domain_error("Read non-integer value for green");
      }
      if (!isValidColor(image[col][row].r))
      {
        ostringstream oss;
        oss << "Invalid color value for green (" << image[col][row].r + ")";
        throw std::domain_error(oss.str());
      }

      ifs >> image[col][row].b;
      if (ifs.fail() && !ifs.eof())
      { // bad input that is not caused by being at the end of the file
        throw std::domain_error("Read non-integer value for blue");
      }
      if (!isValidColor(image[col][row].r))
      {
        ostringstream oss;
        oss << "Invalid color value for blue (" << image[col][row].r + ")";
        throw std::domain_error(oss.str());
      }
    }
  }
  cout << "End loadImage..." << endl;
  return image;
}

void outputImage(string filename, const Pixel *const *image, int width, int height)
{
  cout << "Start outputImage..." << endl;
  // remove code
  // declare/define and open output file stream with filename
  ofstream ofs(filename);
  // ensure file is open
  if (!ofs.is_open())
  {
    throw std::invalid_argument("Error: failed to open output file - " + filename);
  }
  ofs << "P3" << endl;
  ofs << width << " " << height << endl;
  ofs << 255 << endl;
  for (int row = 0; row < height; ++row)
  {
    for (int col = 0; col < width; ++col)
    {
      ofs << image[col][row].r << " ";
      ofs << image[col][row].g << " ";
      ofs << image[col][row].b << " ";
    }
    ofs << endl;
  }
  cout << "End outputImage..." << endl;
}