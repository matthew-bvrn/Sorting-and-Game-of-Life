#include <iostream>
#include <stdlib.h>
#include <exception>
#include <assert.h>
#include <ctime>
#include <chrono>
#include <fstream>
#include <vector>
#include <ios>
#include<cmath>
#include <algorithm>
//#ifndef SORTABLECON_H
//#define SORTABLECON_H

struct IntegerCoordinate
{

	int X, Y;
};

class SortableContainer
{
public:
	virtual	bool cmp(int i, int j) = 0;
	virtual int size() = 0;
	virtual void swap(int i, int j) = 0;
};


class CoordinateArray : public SortableContainer
{
public:
	//constructor
	explicit CoordinateArray(int n) : v(n) {}
	//return values
	int x(int i)
	{

		return v[i].X;
	}
	int y(int i)
	{
		return v[i].Y;
	}
	//access values
	int &operator()(int index, int coord)
	{
		if (coord == 0)
		{

			return v[index].X;
		}
		else
		{

			return v[index].Y;
		}
	}
	//write format
	friend std::ostream &operator<<(std::ostream &os, CoordinateArray v)
	{
		os << "(";
		for (int i = 0; i < v.size() - 1; i++)
		{
			os << "(" << v.x(i) << ", " << v.y(i) << ") ,";
		}
		os << "(" << v.x(v.size() - 1) << ", " << v.y(v.size() - 1) << ")";
		os << ")" << std::endl;
		return os;
	}
	//lexicographic order
	bool cmp(int i, int j)
	{
		if (x(i) < x(j) || (x(i) == x(j) && y(i) < y(j)))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	//Initialise new array of arrays
	CoordinateArray newLife(int X, int Y)
	{
		CoordinateArray newArray(v.size() + 1);
		for (int i = 0; i < newArray.size() - 1; i++)
		{
			newArray(i, 0) = x(i);
			newArray(i, 1) = y(i);
		}
		newArray(newArray.size() - 1, 0) = X;
		newArray(newArray.size() - 1, 1) = Y;
		return newArray;
	}

	//removes element from given coordinate
	CoordinateArray deleter(int i, CoordinateArray current)
	{
		CoordinateArray newArray(current.size() - 1);
		for (int j = 0; j < i; j++)
		{
			newArray(j, 0) = current.x(j);
			newArray(j, 1) = current.y(j);


		}
		for (int j = i; j < current.size() - 1; j++)
		{

			newArray(j, 0) = current.x(j + 1);
			newArray(j, 1) = current.y(j + 1);
		}

		return newArray;
	}


	//swap piecewise
	void swap(int i, int j)
	{
		double swapper = v[i].X;
		v[i].X = v[j].X;
		v[j].X = swapper;
		swapper = v[i].Y;
		v[i].Y = v[j].Y;
		v[j].Y = swapper;
	}

	int size() { return v.size(); }
	//randomise each coordinate separately
	void initialise_random(double xmin, double xmax)
	{
		size_t s = v.size();
		for (size_t i = 0; i < s; i++)
		{
			v[i].X = (xmin + (xmax - xmin)*rand() / static_cast<double>(RAND_MAX));
			v[i].Y = (xmin + (xmax - xmin)*rand() / static_cast<double>(RAND_MAX));
		}
	}

	int inArray(int Xcheck, int Ycheck, int lowerBoundary, int upperBoundary)
	{
		if (upperBoundary - lowerBoundary <= 1)
		{
			if (v[upperBoundary].X == Xcheck && v[upperBoundary].Y == Ycheck)
			{
				return upperBoundary;
			}
			else if (v[lowerBoundary].X == Xcheck && v[lowerBoundary].Y == Ycheck)
			{
				return lowerBoundary;
			}
			else
			{
				return -1;
			}


		}
		int checker = lowerBoundary + ceil((upperBoundary - lowerBoundary) / 2);
		//minor speed improvement to look this up only once.
		int VcheckX = v[checker].X;
		int VcheckY = v[checker].Y;

		if (VcheckX == Xcheck && VcheckY == Ycheck)
		{
			return checker;
		}
		else if (VcheckX < Xcheck || (VcheckX == Xcheck && VcheckY < Ycheck))
		{
			return inArray(Xcheck, Ycheck, checker, upperBoundary);
		}
		else
		{
			return inArray(Xcheck, Ycheck, lowerBoundary, checker);
		}
	}

private:
	std::vector<IntegerCoordinate> v;
};




#ifndef MVECTOR_H
#define MVECTOR_H

#include <vector>



// Class that represents a mathematical vector
class MVector : public SortableContainer
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	//print overload
	friend std::ostream &operator<<(std::ostream &os, MVector v)
	{
		os << "(";
		for (int i = 0; i < v.size() - 1; i++)
		{
			os << v[i] << ", ";
		}
		os << v[v.size() - 1];
		os << ")" << std::endl;
		return os;
	}
	// access element (write) 
	double &operator[](int index)
	{
		if (index > size() - 1)
		{
			
		}
		return v[index];
	}

	//compare size
	bool cmp(int i, int j)
	{
		if (v[i] > v[j])
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	//swap
	void swap(int i, int j)
	{
		double swapper = v[i];
		v[i] = v[j];
		v[j] = swapper;
	}
	void initialise_random(double xmin, double xmax)
	{
		size_t s = v.size();
		for (size_t i = 0; i < s; i++)
		{
			v[i] = (xmin + (xmax - xmin)*rand() / static_cast<double>(RAND_MAX));
		}
	}



	int size() { return v.size(); } // number of elements

private:
	std::vector<double> v;

};

#endif

//print a double
void printi(double i) {
	std::cout << i << std::endl;
}

//print a vector
void printvec(MVector v)
{
	for (int i = 0; i <= v.size() - 1; i++)
	{

		printi(v[i]);
	}
	std::cout << std::endl;
}


//bubblesort implementation
void bubble(SortableContainer &v)
{
	int n = v.size();
	for (int j = 0; j < n; j++) {
	
		for (int i = 0; i < n - j - 1; i++)
		{

			if (v.cmp(i + 1, i)) v.swap(i, i + 1);

		}

	}
}

//confirms whether a vector is sorted
int confirm(MVector v)
{
	for (int i = 0; i < v.size() - 1; i++)
	{
		if (v[i] >= v[i + 1])
		{
			continue;
		}
		else
		{
			return 0;
		}
	}
	return 1;
}


void quick_recursive(SortableContainer &v, int start, int end)
{
	//variable for the index being compared
	int beingChecked = start;

	//index of the pivot
	double comparatorIndex = floor(start + ((end)-start)*rand() / static_cast<double>(RAND_MAX));

	//moves down each time so we know where to insert an element larger than the pivot
	int upperIndex = end;

	//put the pivot at the end
	v.swap(comparatorIndex, end);
	comparatorIndex = end;

	//loop until everything greater than the pivot is above it
	while (beingChecked != comparatorIndex)
	{

		if (v.cmp(comparatorIndex, beingChecked))
		{

			//move pivot down to make space for the element greater than it
			v.swap(comparatorIndex, comparatorIndex - 1);
			comparatorIndex--;

			//if they were next to each other we dont want to swap them back so break
			if (comparatorIndex == beingChecked)
			{
				break;
			}

			//insert said element
			v.swap(beingChecked, comparatorIndex + 1);

		}

		//if it was less than or equal to we can compare the next element
		else
		{
			beingChecked++;
		}

	}


	//recurses until the list is 1 or 0 in length
	if (!(comparatorIndex - 1 - start < 1))
	{

		quick_recursive(v, start, comparatorIndex - 1);
	}

	if (!(end - (comparatorIndex + 1) < 1))
	{

		quick_recursive(v, comparatorIndex + 1, end);
	}


}

//heapsort recursion
void heap_from_root(SortableContainer &v, int i, int n)
{
	int biggest_child;

	//if they have no children, theres nothing to do
	if (2 * i + 1 >= n && 2 * i + 2 >= n)
	{
	}

	//if there's only one child, only need to compare the mother to child
	else if (2 * i + 2 >= n)
	{

		//if child is bigger, swap them
		if (v.cmp(i, 2 * i + 1))
		{
			v.swap(2 * i + 1, i);

		}
	}

	//theres two children
	else
	{
		//find out which child is larger
		if (v.cmp(2 * i + 2, 2 * i + 1))
		{
			biggest_child = 2 * i + 1;
		}
		else
		{
			biggest_child = 2 * i + 2;
		}
		//swap with mother if biggest child is bigger
		if (v.cmp(i, biggest_child))
		{
			v.swap(biggest_child, i);
		}

		heap_from_root(v, biggest_child, n);

	}
}

//heapsort function
void heap(SortableContainer &v)
{
	
	for (int j = v.size() - 1; j >= 0; j--)
	{
		heap_from_root(v, j, v.size());
	}


	for (int n = v.size(); n > 0; n--)
	{
		v.swap(0, n - 1);
		heap_from_root(v, 0, n - 1);





	}
}




void quick(SortableContainer &v) { quick_recursive(v, 0, v.size() - 1); }

void rowInitialise(MVector &v)
{
	v[0] = -2;
	v[1] = -2;
	v[2] = -2;
}

void rowShift(MVector &v)
{
	v[0] = v[1];
	v[1] = v[2];
	v[2] = -2;
}

class Life
{
public:
	CoordinateArray LiveCells;

	explicit Life(int n) : LiveCells(n) {}





	CoordinateArray tick(CoordinateArray changedCells)
	{
		
		//CoordinateArray LiveCellsAndNeighours = LiveCells;
		heap(LiveCells);
		
		//determines the size of the field, so we only iterate over a rectangle that contains neighbours of LiveCells
		int maxX = 0;
		int maxY = 0;
		int minX = 100000;
		int minY = 100000;
		for (int i = 0; i < LiveCells.size(); i++)
		{
			if (LiveCells.x(i) > maxX)
			{
				maxX = LiveCells.x(i);
			}
			if (LiveCells.y(i) > maxY)
			{
				maxY = LiveCells.y(i);
			}
			if (LiveCells.x(i) < minX)
			{
				minX = LiveCells.x(i);
			}
			if (LiveCells.y(i) < minY)
			{
				minY = LiveCells.y(i);
			}


		}
		minX--;
		minY--;
		maxX++;
		maxY++;

		//the updated list so each iteration doesn't modify future iterations. Once the tick is complete this overwrites the array.
		CoordinateArray newLiveCells(0);
		CoordinateArray newchangedCells(0);
		//for every cell in the world
		for (int curY = minY; curY <= maxY; curY++)
		{
			//this is an attempt to optimise the code by reducing the number of inArray lookups
			// we need to do by approx 2/3. Since the script moves from column to column in turn,
			// and the results of inArray never change before the end of the tick,
			//We can store the results in the neighbourhood that the following cell
			// will also have so that it can refer to the saved
			//information rather than running inArray again.

			//Firstly, we initialise each of the data storing vectors with -2, indicating they can be overwritten
			MVector AboveRow(3);
			rowInitialise(AboveRow);
			MVector SameRow(3);
			rowInitialise(SameRow);
			MVector BelowRow(3);
			rowInitialise(BelowRow);

			//for each cell 
			for (int curX = minX; curX <= maxX; curX++)
			{
				
				//check if its alive
				int location = LiveCells.inArray(curX, curY, 0, LiveCells.size() - 1);

				//if it is alive
				if (!(location == -1))
				{

					int aliveNeighbours = 0;

					//count up the squares around it
					for (int yDisplace = -1; yDisplace <= 1; yDisplace++)
					{


						for (int xDisplace = -1; xDisplace <= 1; xDisplace++)
						{
							//in each case, we check whether there is already information saved in the vectors. 
							//If not, we then store the result of inArray. If inArray returns -1 it indicates that
							// there is no cell in the CoordinateArray with the prescribed coordiantes.
							//Therefore we check if it returns something else, and if it does, we add +1 to the aliveNeighbours.
							if (yDisplace == -1)
							{
								if ((BelowRow[xDisplace + 1] == -2))
								{
									BelowRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(BelowRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
							else if (yDisplace == 0)
							{
								if ((SameRow[xDisplace + 1] == -2))
								{
									SameRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(SameRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
							else
							{
								if ((AboveRow[xDisplace + 1] == -2))
								{
									AboveRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(AboveRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
						}


					}


					//if it equals 2 or 3 (subtract 1 since the cell itself is counted), keep it by adding it to the list
					//(I originally copied LiveCells to NewLiveCells and then deleted a cell if the reverse was true,
					//But this proved to be timewasting as it had to repeatedly heapsort NewLiveCells
					if ((aliveNeighbours == 3 || aliveNeighbours == 4))
					{
						newLiveCells = newLiveCells.newLife(curX, curY);

					}
					else
					{
						newchangedCells = newchangedCells.newLife(curX, curY);
					}


				}

				//if it is dead
				else
				{
					int aliveNeighbours = 0;

					//count the living neighbours
					for (int xDisplace = -1; xDisplace <= 1; xDisplace++)
					{
						for (int yDisplace = -1; yDisplace <= 1; yDisplace++)
						{

							if (yDisplace == -1)
							{
								if ((BelowRow[xDisplace + 1] == -2))
								{
									BelowRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(BelowRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
							else if (yDisplace == 0)
							{
								if ((SameRow[xDisplace + 1] == -2))
								{
									SameRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(SameRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
							else
							{
								if ((AboveRow[xDisplace + 1] == -2))
								{
									AboveRow[xDisplace + 1] = LiveCells.inArray(curX + xDisplace,
										curY + yDisplace, 0, LiveCells.size() - 1);
								}
								if (!(AboveRow[xDisplace + 1] == -1))
								{
									aliveNeighbours++;
								}
							}
						}
					}




					//if it has three neighbours, bring it to life
					if (aliveNeighbours == 3)
					{			
						newchangedCells = newchangedCells.newLife(curX, curY);
						newLiveCells = newLiveCells.newLife(curX, curY);
						
					}
				}

				//we then shift the information in these vectors along by one, and free up the space
				//where we dont have information by writing it with a -2.
				//This way the next cell being checked will only need to call inArray 3 times, not 9.
				rowShift(AboveRow);
				rowShift(SameRow);
				rowShift(BelowRow);
			}

		}

		//set LiveCells to the new list

		LiveCells = newLiveCells;
		return newchangedCells;

	}
};

int main()
{
	Life cells(7);


	cells.LiveCells(0, 0) = 1001;
	cells.LiveCells(0, 1) = 1001;

	cells.LiveCells(1, 0) = 1002;
	cells.LiveCells(1, 1) = 1001;

	cells.LiveCells(2, 0) = 1002;
	cells.LiveCells(2, 1) = 1003;

	cells.LiveCells(3, 0) = 1004;
	cells.LiveCells(3, 1) = 1002;

	cells.LiveCells(4, 0) = 1005;
	cells.LiveCells(4, 1) = 1001;

	cells.LiveCells(5, 0) = 1006;
	cells.LiveCells(5, 1) = 1001;

	cells.LiveCells(6, 0) = 1007;
	cells.LiveCells(6, 1) = 1001;


	std::cout << cells.LiveCells;
	CoordinateArray changedCells = cells.LiveCells;
	for (int i = 0; i < 5206; i++)
	{
		changedCells=cells.tick(changedCells);
		std::cout << i << std::endl;

	}




	std::cout << cells.LiveCells;

}
