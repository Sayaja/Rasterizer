#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t; // Time between Update()
float yaw; // Rotation angle
mat3 R; // Rotation matrix
vector<Triangle> triangles;
float focalLength = SCREEN_WIDTH;
vec3 cameraPos( 0, 0, -3.001 );
vec3 currentColor; // Color of the current triangle
vec3 currentNormal; // Normal of current triangle
vec3 currentReflectance; // Reflectance of current triangle
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 14.1f * vec3(1,1,1);
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );
struct Pixel {
int x;
int y;
float zinv;
vec3 pos3d;
};
struct Vertex {
vec3 position;
};


// ----------------------------------------------------------------------------
// FUNCTIONS

void VertexShader(const vec3& v, Pixel& p);
void PixelShader(const Pixel& p);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void InterpolateFloat(float a, float b, vector<float>& result);
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void InterpolateGLM(vec3 a, vec3 b, vector<vec3>& resultGLM);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertexPixels,vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels);
void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);
void Update();
void Draw();

// ----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	// Set values for the rotation matrix
	R[0][0] = cos(yaw); R[0][1] = 0; R[0][2] = sin(yaw);
	R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
	R[2][0] = -sin(yaw); R[2][1] = 0; R[2][2] = cos(yaw);

	// To test ComputePolygonRows
	// vector<ivec2> vertexPixels(3);
	// vertexPixels[0] = ivec2(10, 5);
	// vertexPixels[1] = ivec2( 5,10);
	// vertexPixels[2] = ivec2(15,15);
	// vector<ivec2> leftPixels;
	// vector<ivec2> rightPixels;
	// ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	// for( int row=0; row<leftPixels.size(); ++row ) {
	// 	cout << "Start: ("
	// 	<< leftPixels[row].x << ","
	// 	<< leftPixels[row].y << "). "
	// 	<< "End: ("
	// 	<< rightPixels[row].x << ","
	// 	<< rightPixels[row].y << "). " << endl;
	// }

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void VertexShader(const Vertex& v, Pixel& p) { // Transform vertex from 3d -> 2d
	vec3 P = (v.position - cameraPos) * R;
	p.pos3d = v.position;
	p.zinv = 1.0f / P.z;
	p.x = ((focalLength * P.x) / P.z) + (SCREEN_WIDTH / 2);
	p.y = ((focalLength * P.y) / P.z) + (SCREEN_HEIGHT / 2);
}


// void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) {
// 	int N = result.size();
// 	ivec2 step = ivec2(b-a) / (glm::max(N-1,1));
// 	ivec2 current( a );
// 	for( int i=0; i<N; ++i ) {
// 		result[i] = current;
// 		current += step;
// 	}
// }

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) { // Replaced with this one cuz c++ wouldn't
	float temp;                                               // accept float(glm::max(N-1,1)). Makes it
	int N = result.size();                                    // ~25% slower when rendering
	if (N-1 <= 1) {
		temp = 1;
	} else {
		temp = N - 1;
	}
	float diff0 = b[0] - a[0];
	float step0 = diff0 / temp;
	float diff1 = b[1] - a[1];
	float step1 = diff1 / temp;
	for(int i=0; i<result.size(); ++i) {
		result[i].x = a[0] + i * step0;
		result[i].y = a[1] + i * step1;
	}
}

void InterpolateFloat(float a, float b, vector<float>& result) {
	float diff = b - a; // Get the difference between b and a
	float step = diff / (result.size() - 1); // Calculate the difference between each interpolated value
	for(int i=0; i<result.size(); ++i) {
		result[i] = a + i * step;
	}
	//return result;
}

void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result) {
	float temp;
	int N = result.size();
	if (N-1 <= 1) {
		temp = 1;
	} else {
		temp = N - 1;
	}
	float diff0 = b.x - a.x;
	float step0 = diff0 / temp;
	float diff1 = b.y - a.y;
	float step1 = diff1 / temp;
	float diff2 = b.zinv - a.zinv;
	float step2 = diff2 / temp;
	//vec3 diff3 = b.illumination - a.illumination;
	//vec3 step3 = diff3 / temp;
	float diff3 = b.pos3d[0] - a.pos3d[0];
	float step3 = diff3 / temp;
	float diff4 = b.pos3d[1] - a.pos3d[1];
	float step4 = diff4 / temp;
	float diff5 = b.pos3d[2] - a.pos3d[2];
	float step5 = diff5 / temp;
	for(int i=0; i<result.size(); ++i) {
		result[i].x = a.x + i * step0;
		result[i].y = a.y + i * step1;
		result[i].zinv = a.zinv + i * step2;
		result[i].pos3d[0] = a.pos3d[0] + i * step3;
		result[i].pos3d[1] = a.pos3d[1] + i * step4;
		result[i].pos3d[2] = a.pos3d[2] + i * step5;
	}
}

void InterpolateGLM(vec3 a, vec3 b, vector<vec3>& resultGLM) {
	int temp;
	int N = resultGLM.size();
	if (N-1 <= 1) {
		temp = 1;
	} else {
		temp = N - 1;
	}

	float diff0 = b[0] - a[0];
	float step0 = diff0 / temp;
	float diff1 = b[1] - a[1];
	float step1 = diff1 / temp;
	float diff2 = b[2] - a[2];
	float step2 = diff2 / temp;
	for(int i=0; i<resultGLM.size(); ++i) {
		resultGLM[i].x = a[0] + i * step0;
		resultGLM[i].y = a[1] + i * step1;
		resultGLM[i].z = a[2] + i * step2;
	}
}

// Draw a line between two pixels
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color) {
	ivec2 delta = glm::abs( a - b );
	int pixels = glm::max( delta.x, delta.y ) + 1;
	vector<ivec2> line( pixels );
	Interpolate( a, b, line );
	for (int i=0; i<line.size(); ++i) {
		PutPixelSDL(surface, line[i].x, line[i].y, color);
	}
}

void DrawPolygonEdges(const vector<vec3>& vertices) { // Draw lines between vertices
	// int V = vertices.size();
	// // Transform each vertex from 3D world position to 2D image position:
	// vector<ivec2> projectedVertices( V );
	// for( int i=0; i<V; ++i ) { // Get the projected 2d vertices
	// 	VertexShader( vertices[i], projectedVertices[i] );
	// }
	// // Loop over all vertices and draw the edge from it to the next vertex:
	// for( int i=0; i<V; ++i ) {
	// 	int j = (i+1)%V; // The next vertex
	// 	vec3 color( 1, 1, 1 );
	// 	DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	// }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels) { // Find out which pixels a polygon ocupy
// 1. Find max and min y-value of the polygon
// and compute the number of rows it occupies.
// 2. Resize leftPixels and rightPixels
// so that they have an element for each row.
// 3. Initialize the x-coordinates in leftPixels
// to some really large value and the x-coordinates
// in rightPixels to some really small value.
// 4. Loop through all edges of the polygon and use
// linear interpolation to find the x-coordinate for
// each row it occupies. Update the corresponding
// values in rightPixels and leftPixels.

	int minY = SCREEN_HEIGHT;
	int maxY = 0;
	for (int i=0;i<vertexPixels.size();++i) { // Get the largest and smallest for the polygon
		if (vertexPixels[i].y < minY) {
			minY = vertexPixels[i].y;
		}
		if (maxY < vertexPixels[i].y) {
			maxY = vertexPixels[i].y;
		}
	}
	int numRows = maxY - minY + 1; // Number of rows
	leftPixels.resize(numRows); // The leftmost pixel for each row
	rightPixels.resize(numRows); // The rightmost pixel for each row
	for (int i=0; i<numRows; ++i) {
		leftPixels[i].x = +numeric_limits<int>::max(); // Set starting values for x
		rightPixels[i].x = -numeric_limits<int>::max();
		leftPixels[i].y = minY + i; // Set y-values
		rightPixels[i].y = minY + i;
	}

	for (int i=0;i<vertexPixels.size();++i) { // Loop through every edge of the polygon
		int l = (i+1)%vertexPixels.size(); // The next vertex
		float deltaX = glm::abs( vertexPixels[i].x - vertexPixels[l].x );
		float deltaY = glm::abs( vertexPixels[i].y - vertexPixels[l].y );
		int pixels = glm::max( deltaX, deltaY ) + 1; // Amount of values in the interpolation
		vector<Pixel> line( pixels );
		InterpolatePixel(vertexPixels[i], vertexPixels[l], line);
		for (int j=0;j<line.size();++j) { // For every step in the interpolated line
			for (int k=0;k<numRows;++k) { // For every row in the polygon
				if (int(line[j].y) == leftPixels[k].y) { // Check if they are matching elements
					if (int(line[j].x) < leftPixels[k].x) { // Replace x values
						leftPixels[k].x = int(line[j].x);
						leftPixels[k].zinv = line[j].zinv;
						leftPixels[k].pos3d = line[j].pos3d;
					}
					if (int(line[j].x) > rightPixels[k].x) {
						rightPixels[k].x = int(line[j].x);
						rightPixels[k].zinv = line[j].zinv;
						rightPixels[k].pos3d = line[j].pos3d;
					}
				}
			}
		}

	}
}

void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
	for (int i=0;i<leftPixels.size();++i) {
		vector<float> rowZInv((rightPixels[i].x - leftPixels[i].x) + 1);
		vector<vec3> rowPos3d((rightPixels[i].x - leftPixels[i].x) + 1);

		// leftPixels[i].pos3d.x = 1.0f / leftPixels[i].pos3d.x;
		// leftPixels[i].pos3d.y = 1.0f / leftPixels[i].pos3d.y;;
		// leftPixels[i].pos3d.z = 1.0f / leftPixels[i].pos3d.z;;
		// rightPixels[i].pos3d.x = 1.0f / rightPixels[i].pos3d.x;
		// rightPixels[i].pos3d.y = 1.0f / rightPixels[i].pos3d.y;
		// rightPixels[i].pos3d.z = 1.0f / rightPixels[i].pos3d.z;

		InterpolateFloat(leftPixels[i].zinv, rightPixels[i].zinv, rowZInv); // Only interpolate zinv for faster rendering
		InterpolateGLM(leftPixels[i].pos3d, rightPixels[i].pos3d, rowPos3d); // Interpolate 3d position
		for (int j=leftPixels[i].x; j <= rightPixels[i].x; ++j) {
			if (rowZInv[j-leftPixels[i].x] > depthBuffer[j][leftPixels[i].y]) {
				depthBuffer[j][leftPixels[i].y] = rowZInv[j-leftPixels[i].x];

				//rowPos3d[j-leftPixels[i].x] = rowPos3d[j-leftPixels[i].x] * rowZInv[j-leftPixels[i].x];

				vec3 r = lightPos - rowPos3d[j-leftPixels[i].x];
				float rad = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
				r = r / rad;
				float temp = dot(r, currentNormal);
				if (temp < 0)
				{
					temp = 0;
				}
				vec3 D;
				D[0] = lightPower[0] * (temp / (4 * 3.14*rad*rad));
				D[1] = lightPower[1] * (temp / (4 * 3.14*rad*rad));
				D[2] = lightPower[2] * (temp / (4 * 3.14*rad*rad));
				vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);

				PutPixelSDL(screen, j, leftPixels[i].y, illumination);
			}
		}
	}
}

// Combine all functions to draw the complete polygon
void DrawPolygon(const vector<Vertex>& vertices) {
	int V = vertices.size();
	vector<Pixel> vertexPixels( V );
	for( int i=0; i<V; ++i ) {
		VertexShader( vertices[i], vertexPixels[i] );
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	DrawPolygonRows( leftPixels, rightPixels );
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);

	if( keystate[SDLK_UP] ) {
		// Move camera forward
		vec3 tempVec(0,0,0.2);
		cameraPos = cameraPos + R * tempVec;
	}

	if( keystate[SDLK_DOWN] ) {
		// Move camera backward
		vec3 tempVec(0,0,-0.2);
		cameraPos = cameraPos + R * tempVec;
	}

	if( keystate[SDLK_RIGHT] ) {
		//Rotate camera counterclockwise
		yaw = yaw + float(-(3.14/16));
		R[0][0] = cos(yaw); R[0][1] = 0; R[0][2] = sin(yaw);
		R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
		R[2][0] = -sin(yaw); R[2][1] = 0; R[2][2] = cos(yaw);
	}

	if( keystate[SDLK_LEFT] ) {
		//Rotate camera clockwise
		yaw = yaw + float((3.14/16));
	  R[0][0] = cos(yaw); R[0][1] = 0; R[0][2] = sin(yaw);
		R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
		R[2][0] = -sin(yaw); R[2][1] = 0; R[2][2] = cos(yaw);
	}

	if( keystate[SDLK_w] ) {
		vec3 tempVec(0, 0, 0.2);
		lightPos = lightPos + R * tempVec;
	}

	if( keystate[SDLK_s] ) {
		vec3 tempVec(0, 0, -0.2);
		lightPos = lightPos + R * tempVec;
	}

	if( keystate[SDLK_d] ) {
		vec3 tempVec(0.2, 0, 0);
		lightPos = lightPos + R * tempVec;
	}

	if( keystate[SDLK_a] ) {
		vec3 tempVec(-0.2, 0, 0);
		lightPos = lightPos + R * tempVec;
	}

	if( keystate[SDLK_q] ) {
		vec3 tempVec(0, 0.2, 0);
		lightPos = lightPos + R * tempVec;
	}

	if( keystate[SDLK_e] ) {
		vec3 tempVec(0, -0.2, 0);
		lightPos = lightPos + R * tempVec;
	}
}

void Draw()
{
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	for( int y=0; y<SCREEN_HEIGHT; ++y ) {
		for( int x=0; x<SCREEN_WIDTH; ++x ) {
			depthBuffer[y][x] = 0;
		}
	}

	// Check if DrawLineSDL works
	//ivec2 a (359,257);
	//ivec2 b(101,99);
	//vec3 color2(1,0,0);
	//DrawLineSDL(screen, a, b, color2);

	for( int i=0; i<triangles.size(); ++i )
	{
		currentColor = triangles[i].color; // Get the color of the current triangle
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		//DrawPolygonEdges(vertices);
		DrawPolygon(vertices);

		// for(int v=0; v<3; ++v) {
		// 	ivec2 projPos;
		// 	VertexShader( vertices[v], projPos );
		// 	vec3 color(1,1,1);
		// 	PutPixelSDL( screen, projPos.x, projPos.y, color );
		// }
	}

	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
