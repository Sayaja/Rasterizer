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

const int SCREEN_WIDTH = 400;
const int SCREEN_HEIGHT = 400;
SDL_Surface* screen;
int t; // Time between Update()
float yaw; // Rotation angle
mat3 R; // Camera rotation matrix
mat3 RL; // Light source rotation matrix
vector<Triangle> triangles;
float focalLength = SCREEN_WIDTH;
vec3 cameraPos( 0.0f, 0.0f, -3.001f );
vec3 currentNormal; // Normal of current triangle
vec3 currentReflectance; // Reflectance of current triangle
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; // Inverse depth for each pixel
float lightBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; // Inverse depth from light source for each pixel
vec3 lightPos(1.5f,-0.5f,-4.001f); // Light source position
vec3 lightPower = 100.1f * vec3(1,1,1); // Light source power
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 ); // Set indirect light to constant
struct Pixel { // Information about a pixel
int x;
int y;
int xLight;
int yLight;
float zinv;
float linv;
vec3 pos3d;
};
struct Vertex { // Information about a vertex
vec3 position;
};

// ----------------------------------------------------------------------------
// FUNCTIONS

void VertexShader(const vec3& v, Pixel& p);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void InterpolateFloat(float a, float b, vector<float>& result);
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void InterpolatePixelLight(Pixel a, Pixel b, vector<Pixel>& result);
void InterpolateGLM(vec3 a, vec3 b, vector<vec3>& resultGLM);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels);
void ComputePolygonRowsLight(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels);
void ShadowMapping(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);
void DrawPolygonLight(const vector<Vertex>& vertices);
void Update();
void Draw();

// ----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles ); // Load the scene
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	// Set values for the rotation matrices
	R[0][0] = cos(yaw); R[0][1] = 0; R[0][2] = sin(yaw);
	R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
	R[2][0] = -sin(yaw); R[2][1] = 0; R[2][2] = cos(yaw);
	mat3 RL1;
	mat3 RL2;
	RL1[0][0] = cos(3.14/10); RL1[0][1] = 0; RL1[0][2] = sin(3.14/10);
	RL1[1][0] = 0; RL1[1][1] = 1; RL1[1][2] = 0;
	RL1[2][0] = -sin(3.14/10); RL1[2][1] = 0; RL1[2][2] = cos(3.14/10);
	RL2[0][0] = 1; RL2[0][1] = 0; RL2[0][2] = 0;
	RL2[1][0] = 0; RL2[1][1] = cos(3.14/16); RL2[1][2] = -sin(3.14/16);
	RL2[2][0] = 0; RL2[2][1] = sin(3.14/16); RL2[2][2] = cos(3.14/16);
	RL = RL1 * RL2;

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
	// Transform to camera pixels
	float x = ((focalLength * P.x) / P.z) + (SCREEN_WIDTH / 2);
	float y = ((focalLength * P.y) / P.z) + (SCREEN_HEIGHT / 2);
	if (x > 0) {
		p.x = x + 0.5;
	} else {
		p.x = x - 0.5;
	}
	if (y > 0) {
		p.y = y + 0.5;
	} else {
		p.y = y - 0.5;
	}
	vec3 P2 = (v.position - lightPos) * RL;
	p.linv = 1 / (sqrt(P2[0]*P2[0] + P2[1]*P2[1] + P2[2]*P2[2])); // Inversed
	// Transform to light pixels
	float xLight = ((focalLength * P2.x) / P2.z) + (SCREEN_WIDTH / 2);
	float yLight = ((focalLength * P2.y) / P2.z) + (SCREEN_HEIGHT / 2);
	if (xLight > 0) {
		p.xLight = xLight + 0.5;
	} else {
		p.xLight = xLight - 0.5;
	}
	if (yLight > 0) {
		p.yLight = yLight + 0.5;
	} else {
		p.yLight = yLight - 0.5;
	}
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) { // Interpolate ivec2
	float temp;
	int N = result.size();
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

void InterpolateFloat(float a, float b, vector<float>& result) { // Interpolate float
	float diff = b - a; // Get the difference between b and a
	float step = diff / (result.size() - 1); // Calculate the difference between each interpolated value
	for(int i=0; i<result.size(); ++i) {
		result[i] = a + i * step;
	}
}

void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result) { // Interpolate pixels (excluding light coordinates)
	float temp;
	int N = result.size();
	if (N-1 <= 1) {
		temp = 1;
	} else {
		temp = N - 1;
	}

	a.pos3d = a.pos3d * a.zinv; // Perspective correct interpolation
	b.pos3d = b.pos3d * b.zinv;

	float diff0 = b.x - a.x;
	float step0 = diff0 / temp;
	float diff1 = b.y - a.y;
	float step1 = diff1 / temp;
	float diff2 = b.zinv - a.zinv;
	float step2 = diff2 / temp;
	float diff3 = b.pos3d[0] - a.pos3d[0];
	float step3 = diff3 / temp;
	float diff4 = b.pos3d[1] - a.pos3d[1];
	float step4 = diff4 / temp;
	float diff5 = b.pos3d[2] - a.pos3d[2];
	float step5 = diff5 / temp;
	for(int i=0; i<result.size(); ++i) {
		float x = a.x + i * step0;
		float y = a.y + i * step1;
		if (x > 0) {
			result[i].x = x + 0.5;
		} else {
			result[i].x = x - 0.5;
		}
		if (y > 0) {
			result[i].y = y + 0.5;
		} else {
			result[i].y = y - 0.5;
		}
		result[i].zinv = a.zinv + i * step2;
		result[i].pos3d[0] = (a.pos3d[0] + i * step3) / result[i].zinv;
		result[i].pos3d[1] = (a.pos3d[1] + i * step4) / result[i].zinv;
		result[i].pos3d[2] = (a.pos3d[2] + i * step5) / result[i].zinv;
	}
}

void InterpolatePixelLight(Pixel a, Pixel b, vector<Pixel>& result) { // Interpolate (excluding camera coordinates)
	float temp;
	int N = result.size();
	if (N-1 <= 1) {
		temp = 1;
	} else {
		temp = N - 1;
	}

	a.pos3d = a.pos3d * a.linv; // Perspective correct interpolation
	b.pos3d = b.pos3d * b.linv;

	float diff0 = b.xLight - a.xLight;
	float step0 = diff0 / temp;
	float diff1 = b.yLight - a.yLight;
	float step1 = diff1 / temp;
	float diff2 = b.linv - a.linv;
	float step2 = diff2 / temp;
	float diff3 = b.pos3d[0] - a.pos3d[0];
	float step3 = diff3 / temp;
	float diff4 = b.pos3d[1] - a.pos3d[1];
	float step4 = diff4 / temp;
	float diff5 = b.pos3d[2] - a.pos3d[2];
	float step5 = diff5 / temp;
	for(int i=0; i<result.size(); ++i) {
		float xLight = a.xLight + i * step0;
		float yLight = a.yLight + i * step1;
		if (xLight > 0) {
			result[i].xLight = xLight + 0.5;
		} else {
			result[i].xLight = xLight - 0.5;
		}
		if (yLight > 0) {
			result[i].yLight = yLight + 0.5;
		} else {
			result[i].yLight = yLight - 0.5;
		}
		result[i].linv = a.linv + i * step2;
		result[i].pos3d[0] = (a.pos3d[0] + i * step3) / result[i].linv;
		result[i].pos3d[1] = (a.pos3d[1] + i * step4) / result[i].linv;
		result[i].pos3d[2] = (a.pos3d[2] + i * step5) / result[i].linv;
	}
}

void InterpolateGLM(vec3 a, vec3 b, vector<vec3>& resultGLM) { // Interpolate vec3
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

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels) { // Find out which pixels a polygon ocupy frin the camera's pow
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
		int deltaX = glm::abs( vertexPixels[i].x - vertexPixels[l].x );
		int deltaY = glm::abs( vertexPixels[i].y - vertexPixels[l].y );
		int pixels = glm::max( deltaX, deltaY ) + 1; // Amount of values in the interpolation
		vector<Pixel> line( pixels );
		InterpolatePixel(vertexPixels[i], vertexPixels[l], line);
		for (int j=0;j<line.size();++j) { // For every step in the interpolated line
			for (int k=0;k<numRows;++k) { // For every row in the polygon
				if (int(line[j].y) == leftPixels[k].y) { // Check if they are matching elements
					if (int(line[j].x) < leftPixels[k].x) { // Replace x values
						leftPixels[k].x = int(line[j].x);
						leftPixels[k].zinv = line[j].zinv; // Set zinv and pos3d aswell
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

void ComputePolygonRowsLight(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
vector<Pixel>& rightPixels) { // Find out which pixels a polygon ocupy from the light's pow

	int minY = SCREEN_HEIGHT;
	int maxY = 0;
	for (int i=0;i<vertexPixels.size();++i) { // Get the largest and smallest for the polygon
		if (vertexPixels[i].yLight < minY) {
			minY = vertexPixels[i].yLight;
		}
		if (maxY < vertexPixels[i].yLight) {
			maxY = vertexPixels[i].yLight;
		}
	}
	int numRows = maxY - minY + 1; // Number of rows
	leftPixels.resize(numRows); // The leftmost pixel for each row
	rightPixels.resize(numRows); // The rightmost pixel for each row
	for (int i=0; i<numRows; ++i) {
		leftPixels[i].xLight = +numeric_limits<int>::max(); // Set starting values for x
		rightPixels[i].xLight = -numeric_limits<int>::max();
		leftPixels[i].yLight = minY + i; // Set y-values
		rightPixels[i].yLight = minY + i;
	}

	for (int i=0;i<vertexPixels.size();++i) { // Loop through every edge of the polygon
		int l = (i+1)%vertexPixels.size(); // The next vertex
		int deltaX = glm::abs( vertexPixels[i].xLight - vertexPixels[l].xLight );
		int deltaY = glm::abs( vertexPixels[i].yLight - vertexPixels[l].yLight );
		int pixels = glm::max( deltaX, deltaY ) + 1; // Amount of values in the interpolation
		vector<Pixel> line( pixels );
		InterpolatePixelLight(vertexPixels[i], vertexPixels[l], line);
		for (int j=0;j<line.size();++j) { // For every step in the interpolated line
			for (int k=0;k<numRows;++k) { // For every row in the polygon
				if (int(line[j].yLight) == leftPixels[k].yLight) { // Check if they are matching elements
					if (int(line[j].xLight) < leftPixels[k].xLight) { // Replace x values
						leftPixels[k].xLight = int(line[j].xLight);
						leftPixels[k].linv = line[j].linv; // Set linv and pos3d aswell
						leftPixels[k].pos3d = line[j].pos3d;
					}
					if (int(line[j].xLight) > rightPixels[k].xLight) {
						rightPixels[k].xLight = int(line[j].xLight);
						rightPixels[k].linv = line[j].linv;
						rightPixels[k].pos3d = line[j].pos3d;
					}
				}
			}
		}
	}
}

// Fill the lightBuffer matrix
void ShadowMapping(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
	for (int i=0;i<leftPixels.size();++i) { // For every row
		vector<Pixel> row((rightPixels[i].xLight - leftPixels[i].xLight) + 1);

		InterpolatePixelLight(leftPixels[i], rightPixels[i], row); // Interpolate over the row
		for (int j=leftPixels[i].xLight; j <= rightPixels[i].xLight; ++j) { // For every pixel in the row
			if (row[j-leftPixels[i].xLight].linv >= lightBuffer[j][leftPixels[i].yLight]) { // Check if the pixel is closer
				lightBuffer[j][leftPixels[i].yLight] = row[j-leftPixels[i].xLight].linv;
				// Comment in this to see the scene from the light's pow
				//PutPixelSDL(screen, j, leftPixelsLight[i].yLight, currentReflectance);
				}
		}
	}
}

void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
	for (int i=0;i<leftPixels.size();++i) { // For every row
		vector<Pixel> row((rightPixels[i].x - leftPixels[i].x) + 1);

		InterpolatePixel(leftPixels[i], rightPixels[i], row); // Interpolate over the row
		for (int j=leftPixels[i].x; j <= rightPixels[i].x; ++j) { // For every pixel in the row
			if (row[j-leftPixels[i].x].zinv >= depthBuffer[j][leftPixels[i].y]) { // Check if the pixel is closer
				depthBuffer[j][leftPixels[i].y] = row[j-leftPixels[i].x].zinv;

				vec3 P = (row[j-leftPixels[i].x].pos3d - lightPos) * RL;
				int xLight = ((focalLength * P.x) / P.z) + (SCREEN_WIDTH / 2) + 0.5;
				int yLight = ((focalLength * P.y) / P.z) + (SCREEN_HEIGHT / 2) + 0.5;
				float lightDepth = 1.0f / (sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2])); // Inverse distance

				// 0.005
				if ((lightDepth + 0.003) >= lightBuffer[xLight][yLight]) { // Check if it's visible from the light
					// Calculate the illumination
					vec3 r = lightPos - row[j-leftPixels[i].x].pos3d;
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
				else {
					vec3 illumination = currentReflectance * indirectLightPowerPerArea;
					PutPixelSDL(screen, j, leftPixels[i].y, illumination);
				}
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
	ComputePolygonRows( vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows( leftPixels, rightPixels );
}

// Draw the shadow map
void DrawPolygonLight(const vector<Vertex>& vertices) {
	int V = vertices.size();
	vector<Pixel> vertexPixels( V );
	for( int i=0; i<V; ++i ) {
		VertexShader( vertices[i], vertexPixels[i] );
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRowsLight( vertexPixels, leftPixels, rightPixels);
	ShadowMapping(leftPixels, rightPixels);
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

	if (keystate[SDLK_z]) {
		//move camera left
		vec3 temp(-0.2, 0, 0);
		cameraPos = cameraPos + R * temp;
	}

	if (keystate[SDLK_x]) {
		//move camera right
		vec3 temp(0.2, 0, 0);
		cameraPos = cameraPos + R * temp;
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
			lightBuffer[y][x] = 0;
		}
	}

	// First draw the shadow map and fill the lightBuffer matrix
	for( int i=0; i<triangles.size(); ++i )
	{
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		DrawPolygonLight(vertices);
	}

	// Then draw the scene from the camera's pow
	for( int i=0; i<triangles.size(); ++i )
	{
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		DrawPolygon(vertices);
	}

	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
