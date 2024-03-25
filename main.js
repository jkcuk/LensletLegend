
// This code is based on three.js, which comes with the following license:
//
// The MIT License
//
// Copyright © 2010-2024 three.js authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
import * as THREE from 'three';

import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let name = 'LensletLegend';
let scene;
let aspectRatioVideoFeedU = 4.0/3.0;
let aspectRatioVideoFeedE = 4.0/3.0;
let renderer;
let videoFeedU, videoFeedE;	// feeds from user/environment-facing cameras
let camera;
let controls;
let raytracingSphere;
let raytracingSphereShaderMaterial;
	
// Nokia HR20, according to https://www.camerafv5.com/devices/manufacturers/hmd_global/nokia_xr20_ttg_0/
let fovVideoFeedU = 67.3;	// (user-facing) camera
let fovVideoFeedE = 68.3;	// (environment-facing) camera
let fovScreen = 68;

let cameraLensDistance = 10.0;
let raytracingSphereRadius = 20.0;
let offsetFromConfocal = 0.0;
let deltaPeriod = 0.0;

// the status text area
let status = document.createElement('div');
let statusTime;	// the time the last status was posted

// the info text area
let info = document.createElement('div');
let infoTime;	// the time the last info was posted

let gui;

// true if stored photo is showing
let showingStoredPhoto = false;
let storedPhoto;
let storedPhotoDescription;
let storedPhotoInfoString;

// from https://github.com/4nt0nio/jpegcam
const click = new Audio('./shutter.mp3');

// uncomment the foolowing lines, and 
// stats.begin();
// stats.end();
// in animate(), to show fps stats
// import Stats from 'stats.js'
// var stats = new Stats();
// stats.showPanel( 0 ); // 0: fps, 1: ms, 2: mb, 3+: custom
// document.body.appendChild( stats.dom );

init();
animate();

function init() {
	// create the info element first so that any problems can be communicated
	createStatus();

	scene = new THREE.Scene();
	// scene.background = new THREE.Color( 'skyblue' );
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera = new THREE.PerspectiveCamera( fovScreen, windowAspectRatio, 0.1, 2*raytracingSphereRadius + 1 );
	camera.position.z = cameraLensDistance;
	screenChanged();
	
	renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );
	// document.getElementById('livePhoto').appendChild( renderer.domElement );

	createVideoFeeds();

	addRaytracingSphere();

	// user interface

	// handle device orientation
	// window.addEventListener("deviceorientation", handleOrientation, true);

	// handle window resize
	window.addEventListener("resize", onWindowResize, false);

	// share button functionality
	document.getElementById('takePhotoButton').addEventListener('click', takePhoto);

	// toggle fullscreen button functionality
	document.getElementById('fullscreenButton').addEventListener('click', toggleFullscreen);

	// back button functionality
	document.getElementById('backButton').addEventListener('click', showLivePhoto);
	document.getElementById('backButton').style.visibility = "hidden";

	// share button
	document.getElementById('shareButton').addEventListener('click', share);
	document.getElementById('shareButton').style.visibility = "hidden";

	// delete button
	document.getElementById('deleteButton').addEventListener('click', deleteStoredPhoto);
	document.getElementById('deleteButton').style.visibility = "hidden";

	// hide the thumbnail for the moment
	document.getElementById('storedPhotoThumbnail').addEventListener('click', showStoredPhoto);
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('storedPhoto').addEventListener('click', showLivePhoto);
	document.getElementById('storedPhoto').style.visibility = "hidden";
	showingStoredPhoto = false;

	// handle screen-orientation (landscape/portrait) change
	screen.orientation.addEventListener( "change", recreateVideoFeeds );

	addOrbitControls();

	// the controls menu
	createGUI();

	createInfo();
	refreshInfo();
}

function animate() {
	requestAnimationFrame( animate );

	// stats.begin();

	if(!showingStoredPhoto) {
		// update uniforms
		updateUniforms();

		renderer.render( scene,  camera );
	}

	// stats.end();
}

function updateUniforms() {
	// the tangents for the environment-facing camera video feed
	let tanHalfFovHE, tanHalfFovVE;
	if(aspectRatioVideoFeedE > 1.0) {
		// horizontal orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)/aspectRatioVideoFeedE;
	} else {
		// vertical orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)*aspectRatioVideoFeedE;
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.tanHalfFovHE.value = tanHalfFovHE;
	raytracingSphereShaderMaterial.uniforms.tanHalfFovVE.value = tanHalfFovVE;

	// the tangents for the user-facing camera video feed
	let tanHalfFovHU, tanHalfFovVU;
	if(aspectRatioVideoFeedU > 1.0) {
		// horizontal orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)/aspectRatioVideoFeedU;
	} else {
		// vertical orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)*aspectRatioVideoFeedU;
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.tanHalfFovHU.value = tanHalfFovHU;
	raytracingSphereShaderMaterial.uniforms.tanHalfFovVU.value = tanHalfFovVU;

	// calculate the separation between the two arrays, s = f1 + f2 + offsetFromConfocal
	let s = 
		raytracingSphereShaderMaterial.uniforms.focalLength1.value + 
		raytracingSphereShaderMaterial.uniforms.focalLength2.value +
		offsetFromConfocal;
	// arrange them symmetrically around z=0
	raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z = 0;	// +0.5*s;
	raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z = -s;	// -0.5*s;	// - 0.0001;

	// set the array periods
	raytracingSphereShaderMaterial.uniforms.period2.value = raytracingSphereShaderMaterial.uniforms.period1.value + deltaPeriod;
}

function createVideoFeeds() {
	// create the video stream for the user-facing camera first, as some devices (such as my iPad), which have both cameras,
	// but can (for whatever reason) only have a video feed from one at a time, seem to go with the video stream that was
	// created last, and as the standard view is looking "forward" it is preferable to see the environment-facing camera.
	videoFeedU = document.getElementById( 'videoFeedU' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// user-facing camera
		const constraintsU = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'user'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsU ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedU.srcObject = stream;
			videoFeedU.play();

			videoFeedU.addEventListener("playing", () => {
				aspectRatioVideoFeedU = videoFeedU.videoWidth / videoFeedU.videoHeight;
				updateUniforms();
				postStatus(`User-facing(?) camera resolution ${videoFeedU.videoWidth} &times; ${videoFeedU.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access user-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}

	videoFeedE = document.getElementById( 'videoFeedE' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// environment-facing camera
		const constraintsE = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'environment'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsE ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedE.srcObject = stream;
			videoFeedE.play();

			videoFeedE.addEventListener("playing", () => {
				aspectRatioVideoFeedE = videoFeedE.videoWidth / videoFeedE.videoHeight;
				updateUniforms();
				postStatus(`Environment-facing(?) camera resolution ${videoFeedE.videoWidth} &times; ${videoFeedE.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access environment-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}
}

/** create raytracing phere */
function addRaytracingSphere() {
	const videoFeedUTexture = new THREE.VideoTexture( videoFeedU );
	const videoFeedETexture = new THREE.VideoTexture( videoFeedE );
	videoFeedUTexture.colorSpace = THREE.SRGBColorSpace;
	videoFeedETexture.colorSpace = THREE.SRGBColorSpace;

	// the sphere surrouning the camera in all directions
	const geometry = 
		new THREE.SphereGeometry( raytracingSphereRadius );
	raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
		side: THREE.DoubleSide,
		// wireframe: true,
		uniforms: { 
			visible1: { value: true },
			period1: { value: 0.4 },
			alpha1: { value: 0.0 },
			focalLength1: { value: 1.0 },
			centreOfArray1: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			radius1: { value: 5.0 },	// radius of array 1
			visible2: { value: true },
			period2: { value: 0.4 },
			alpha2: { value: 0 },
			focalLength2: { value: -0.5 },
			centreOfArray2: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			radius2: { value: 5.0 },	// radius of array 2
			videoFeedUTexture: { value: videoFeedUTexture }, 
			videoFeedETexture: { value: videoFeedETexture }, 
			// cameraLensDistance: { value: cameraLensDistance },
			tanHalfFovHU: { value: 1.0 },
			tanHalfFovVU: { value: 1.0 },
			tanHalfFovHE: { value: 1.0 },
			tanHalfFovVE: { value: 1.0 }
		},
		vertexShader: `
			varying vec3 intersectionPoint;

			void main()	{
				// projectionMatrix, modelViewMatrix, position -> passed in from Three.js
				intersectionPoint = position.xyz;
  				gl_Position = projectionMatrix
					* modelViewMatrix
					* vec4(position, 1.0);
			}
		`,
		fragmentShader: `
			precision highp float;

			varying vec3 intersectionPoint;
			
			// lenslet array 1
			uniform bool visible1;
			uniform float alpha1;	// rotation angle of array 1
			uniform float period1;	// period of array 1
			uniform float focalLength1;	// focal length of array 1
			uniform vec3 centreOfArray1;	// centre of array 1,  and principal point of lenslet (0, 0)
			uniform float radius1;	// radius of array 1

			// lenslet array 2
			uniform bool visible2;
			uniform float alpha2;	// rotation angle of array 2
			uniform float period2;	// period of array 2
			uniform float focalLength2;	// focal length of array 2
			uniform vec3 centreOfArray2;	// centre of array 2,  and principal point of lenslet (0, 0)
			uniform float radius2;	// radius of array 2

			// video feed from user-facing camera
			uniform sampler2D videoFeedUTexture;
			uniform float tanHalfFovHU;
			uniform float tanHalfFovVU;

			// video feed from environment-facing camera
			uniform sampler2D videoFeedETexture;
			uniform float tanHalfFovHE;
			uniform float tanHalfFovVE;

			// rotate the 2D vector v by the angle alpha (in radians)
			// from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
			vec2 rotate(vec2 v, float alpha) {
				float s = sin(alpha);
				float c = cos(alpha);
				mat2 m = mat2(c, s, -s, c);
				return m * v;
			}

			// Calculate the light-ray direction after transmission through a lens.
			// d is a 2D vector containing the transverse components of the incident light-ray direction,
			// "normalised" such that the longitudinal component of d is 1;
			// r is a 2D vector containing the transverse components of the vector from the principal point
			// to the intersection point;
			// f is the focal length;
			// returns a 2D vector containing the transverse components of the outgoing light-ray direction,
			// "normalised" such that the longitudinal component is 1.
			vec2 lensDeflect(vec2 d, vec2 r, float f) {
				// dPrime propto d/d_z - r/f, where r = I-P, and I = intersection point, P = principal point
				return d - r/f;
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lens.
			// The (ideal thin) lens, of focal length f, is in a z plane through centreOfLens.
			// It is circular, with the given radius, centred on centreOfLenss.
			void passThroughLens(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfLens, 
				float radius,
				float focalLength
			) {
				// "normalised" version of d, scaled such that the z component is 1
				vec3 d1 = d/d.z;

				// calculate the intersection point with the lens
				float deltaZ = centreOfLens.z - p.z;
				p = p + d1*deltaZ;

				if(deltaZ*d.z < 0.0) {
					// the ray actually has to travel *backwards* -- make the lens red
					b *= vec4(0.7, 0.3, 0.3, 1);
				}

				// does the intersection point lie with in the radius?
				vec2 r = p.xy - centreOfLens.xy;
				float r2 = dot(r, r);	// length squared of vector r
				if(r2 < radius*radius) {
					// the intersection point lies inside the radius, so the lens does something to the ray

					// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
					d = vec3(lensDeflect(-d1.xy, r, focalLength), sign(d.z));

					// lower the brightness factor, giving the light a blue tinge
					b *= vec4(0.9, 0.9, 0.99, 1);
				} 
			}

			float findLensletCentreCoordinate(float u, float uPeriod) {
				return uPeriod*floor(u/uPeriod+0.5);
			}

			// Find the principal point of the nearest lenslet in a rectangular lenslet array.
			// r is a 2D vector containing the transverse components of the vector from pp00 to the
			// intersection point;
			// pp00 is the principal point of the lenslet with indices (0, 0) -- sort of the array centre;
			// alpha is the angle (in radians) by which the array is rotated
			// period1 and period2 are the periods in the (rotated) x and y directions
			vec2 findNearestPrincipalPoint(vec2 r, vec2 pp00, float alpha, float period1, float period2) {
				vec2 rr = rotate(r, -alpha);
				vec2 lensletCentreST = vec2(
					findLensletCentreCoordinate(rr.x, period1),
					findLensletCentreCoordinate(rr.y, period2)
				);
				return rotate(lensletCentreST, alpha) + pp00;
			}

			vec2 lensletArrayDeflect(vec2 d, vec2 intersectionPoint, vec2 principalPoint00, float alpha, float period, float focalLength) {
				vec2 r = intersectionPoint - principalPoint00;
				vec2 principalPoint = findNearestPrincipalPoint(
					intersectionPoint - principalPoint00, 
					principalPoint00, alpha, period, period
				);
				// light-ray direction after array
				return lensDeflect(d, intersectionPoint - principalPoint, focalLength);
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lenslet array.
			// The lenslet array is in a z plane through centreOfArray; it is simulated as ideal thin lenses.
			// The lenslets, all of focal length f, are arranged in a square array of the given period, 
			// rotated by alpha around the z axis.  The component is circular, with the given radius, centred on centreOfArray.
			void passThroughLensletArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfArray, 
				float radius,
				float alpha,
				float period,
				float focalLength
			) {
				// transverse part of d, "normalised" such that abs(d.z) = 1
				vec2 d2 = d.xy / abs(d.z);

				// "normalised" version of d, scaled such that the z component is +1
				vec3 d1 = d/d.z;

				// calculate the intersection point with that array
				float deltaZ = centreOfArray.z - p.z;
				p = p + d1*deltaZ;

				// does the intersection point lie with in the radius?
				vec2 r = p.xy - centreOfArray.xy;
				float r2 = dot(r, r);	// length squared of vector r
				if(r2 < radius*radius) {
					// the intersection point lies inside the radius, so the component does something to the ray

					// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
					d = vec3(lensletArrayDeflect(d2, p.xy, centreOfArray.xy, alpha, period, focalLength), sign(d.z));

					// // which direction is the ray passing through it?
					// if(d.z < 0.0) {
					// 	// the ray is passing through array 1 in the direction for which it is designed; deal with it accordingly
					// lower the brightness factor, giving the light a blue tinge
					b *= vec4(0.9, 0.9, 0.99, 1);
					// } else {
					// 	// the ray is passing through the array in the opposite direction for which it is designed; make it red
					// 	b *= vec4(0.7, 0.3, 0.3, 1);
					// }
				
					if(deltaZ*d.z < 0.0) {
						// the ray actually has to travel *backwards* -- make the array red
						b *= vec4(0.7, 0.3, 0.3, 1);
					}
				} 
			}

			void main() {
				// first calculate the current light-ray direction:
				// the camera pinhole is positioned at cameraPosition,
				// the intersection point with the sphere is intersectionPoint,
				// so the "backwards" ray direction from the camera to the intersection point is
				//   d = intersectionPoint - cameraPosition
				vec3 d = intersectionPoint - cameraPosition;

				// the current ray start position; start at the camera
				vec3 p = cameraPosition;

				// current brightness factor; this will multiply the colour at the end
				vec4 b = vec4(1.0, 1.0, 1.0, 1.0);

				if(d.z < 0.0) {
					// the ray is travelling "forwards", in the (-z) direction
					if(visible1) passThroughLensletArray(p, d, b, centreOfArray1,  radius1, alpha1, period1, focalLength1);
					if(visible2) passThroughLensletArray(p, d, b, centreOfArray2,  radius2, alpha2, period2, focalLength2);
					// if(visible1) passThroughLens(p, d, b, centreOfArray1,  radius1, focalLength1);
					// if(visible2) passThroughLens(p, d, b, centreOfArray2,  radius2, focalLength2);
				} else {
					// the ray is travelling "backwards", in the (+z) direction
					if(visible2) passThroughLensletArray(p, d, b, centreOfArray2,  radius2, alpha2, period2, focalLength2);
					if(visible1) passThroughLensletArray(p, d, b, centreOfArray1,  radius1, alpha1, period1, focalLength1);
				}

				// does the ray intersect the (infinitely distant) camera image whose angular width and height is
				// given by arctan(2*tanHalfFovH) and arctan(2*tanHalfFovV?
				vec3 d1 = d/abs(d.z);
				if(d.z < 0.0) {
					// forwards-facing ray
					if((abs(d1.x) < tanHalfFovHE) && (abs(d1.y) < tanHalfFovVE))
						// yes, the ray intersects the image; take the pixel colour from the camera's video feed
						gl_FragColor = texture2D(videoFeedETexture, vec2(0.5+0.5*d1.x/tanHalfFovHE, 0.5+0.5*d1.y/tanHalfFovVE));
					else gl_FragColor = vec4(1, 1, 1, 1.0);	// white
				} else {
					// backwards-facing ray
					if((abs(d1.x) < tanHalfFovHU) && (abs(d1.y) < tanHalfFovVU))
						// yes, the ray intersects the image; take the pixel colour from the camera's video feed
						gl_FragColor = texture2D(videoFeedUTexture, vec2(0.5-0.5*d1.x/tanHalfFovHU, 0.5+0.5*d1.y/tanHalfFovVU));
					else gl_FragColor = vec4(1, 0, 0, 1.0);	// red
				}

				// finally, multiply by the brightness factor
				gl_FragColor *= b;
			}
		`
	});
	raytracingSphere = new THREE.Mesh( geometry, raytracingSphereShaderMaterial ); 
	scene.add( raytracingSphere );
}

// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_animation_skinning_additive_blending.html
function createGUI() {
	// const 
	gui = new GUI();
	// gui.hide();

	const params1 = {
		'Visible': raytracingSphereShaderMaterial.uniforms.visible1.value,
		'Focal length, <i>f</i><sub>1</sub>': raytracingSphereShaderMaterial.uniforms.focalLength1.value,
		'Period, <i>p</i><sub>1</sub>': raytracingSphereShaderMaterial.uniforms.period1.value,
		'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha1.value / Math.PI * 180.
	};
	const params2 = {
		'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
		'Focal length, <i>f</i><sub>2</sub>': raytracingSphereShaderMaterial.uniforms.focalLength2.value,
		'&Delta;<sub>period</sub>, <i>p</i><sub>2</sub> - <i>p</i><sub>1</sub>': deltaPeriod,
		'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
		'Offset from confocal': offsetFromConfocal
	};
	const params = {
		// 'Swap arrays': swapArrays,
		'screen (&deg;)': fovScreen,
		'env.-facing cam. (&deg;)': fovVideoFeedE,
		'user-facing cam. (&deg;)': fovVideoFeedU,
		'Point (virtual) cam. forward (in -<b>z</b> direction)': pointForward,
		'Toggle info visibility': toggleInfoVisibility,
		'Restart video streams': function() { 
			recreateVideoFeed(); 
			postStatus("Restarting video stream");
		}
	}

	const folderArray1 = gui.addFolder( 'Lenslet array 1 (near; in plane <i>z</i><sub>1</sub> = 0)' );

	folderArray1.add( params1, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible1.value = v; } );
	folderArray1.add( params1, 'Focal length, <i>f</i><sub>1</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.focalLength1.value = f; } );
	folderArray1.add( params1, 'Period, <i>p</i><sub>1</sub>', 0.01, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period1.value = p; } );
	folderArray1.add( params1, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha1.value = alpha/180.0*Math.PI; } );

	const folderArray2 = gui.addFolder( 'Lenslet array 2 (far)' );

	folderArray2.add( params1, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible2.value = v; } );
	folderArray2.add( params2, 'Focal length, <i>f</i><sub>2</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.focalLength2.value = f; } );
	folderArray2.add( params2, '&Delta;<sub>period</sub>, <i>p</i><sub>2</sub> - <i>p</i><sub>1</sub>', -0.1, 0.1).onChange( (p) => { deltaPeriod = p; } );
	folderArray2.add( params2, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha2.value = alpha/180.0*Math.PI; } );
	folderArray2.add( params2, 'Offset from confocal', -0.25, 0.25).onChange( (o) => { offsetFromConfocal = o; } );

	const folderFOV = gui.addFolder( 'Fields of view' );
	folderFOV.add( params, 'screen (&deg;)', 10, 170, 1).onChange( setScreenFOV );
	folderFOV.add( params, 'env.-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedE = fov; updateUniforms(); });   
	folderFOV.add( params, 'user-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedU = fov; updateUniforms(); });   
	folderFOV.close();

	const folderSettings = gui.addFolder( 'Other controls' );
	folderSettings.add( params, 'Point (virtual) cam. forward (in -<b>z</b> direction)');
	folderSettings.add( params, 'Toggle info visibility');
	folderSettings.add( params, 'Restart video streams');
	folderSettings.close();
}

/**
 * @param {*} fov	The larger of the camera's horizontal and vertical FOV, in degrees
 * 
 * Set the larger FOV of the screen/window to fov.
 * 
 * Depending on the screen/window's FOV, fov is either the horizontal fov (if screen width > screen height)
 * or the vertical fov (if screen width < screen height).
 */
function setScreenFOV(fov) {
	fovScreen = fov;

	screenChanged();
}

// function swapArrays() {
// 	const visible3 = lensletArrayShaderMaterial.uniforms.visible1.value;
// 	const focalLength3 = lensletArrayShaderMaterial.uniforms.focalLength1.value;
// 	const period3 = lensletArrayShaderMaterial.uniforms.period1.value;
// 	const alpha3 = lensletArrayShaderMaterial.uniforms.alpha1.value;

// 	lensletArrayShaderMaterial.uniforms.visible1.value = lensletArrayShaderMaterial.uniforms.visible2.value;
// 	lensletArrayShaderMaterial.uniforms.focalLength1.value = lensletArrayShaderMaterial.uniforms.focalLength2.value;
// 	lensletArrayShaderMaterial.uniforms.period1.value = lensletArrayShaderMaterial.uniforms.period2.value;
// 	lensletArrayShaderMaterial.uniforms.alpha1.value = lensletArrayShaderMaterial.uniforms.alpha2.value;

// 	lensletArrayShaderMaterial.uniforms.visible2.value = visible3;
// 	lensletArrayShaderMaterial.uniforms.focalLength2.value = focalLength3;
// 	lensletArrayShaderMaterial.uniforms.period2.value = period3;
// 	lensletArrayShaderMaterial.uniforms.alpha2.value = alpha3;
// }

/** 
 * Reset the aspect ratio and FOV of the virtual cameras.
 * 
 * Call if the window size has changed (which also happens when the screen orientation changes)
 * or if camera's FOV has changed
 */
function screenChanged() {
	// alert(`new window size ${window.innerWidth} x ${window.innerHeight}`);

	// in case the screen size has changed
	if(renderer) renderer.setSize(window.innerWidth, window.innerHeight);

	// if the screen orientation changes, width and height swap places, so the aspect ratio changes
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera.aspect = windowAspectRatio;

	// fovS is the screen's horizontal or vertical FOV, whichever is greater;
	// re-calculate the camera FOV, which is the *vertical* fov
	let verticalFOV;
	if(windowAspectRatio > 1.0) {
		// fovS is horizontal FOV; convert to get correct vertical FOV
		verticalFOV = 2.0*Math.atan(Math.tan(0.5*fovScreen*Math.PI/180.0)/windowAspectRatio)*180.0/Math.PI;
	} else {
		// fovS is already vertical FOV
		verticalFOV = fovScreen;
	}
	camera.fov = verticalFOV;

	// make sure the camera changes take effect
	camera.updateProjectionMatrix();
}

function  pointForward() {
	let r = camera.position.length();
	camera.position.x = 0;
	camera.position.y = 0;
	camera.position.z = r;
	controls.update();
	postStatus('Pointing camera forwards (in -<b>z</b> direction)');
}

function onWindowResize() {
	screenChanged();
	postStatus(`window size ${window.innerWidth} &times; ${window.innerHeight}`);	// debug
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeeds() {
	// stop current video streams...
	videoFeedE.srcObject.getTracks().forEach(function(track) { track.stop(); });
	videoFeedU.srcObject.getTracks().forEach(function(track) { track.stop(); });

	// ... and re-create new ones, hopefully of the appropriate size
	createVideoFeeds();
}

function addOrbitControls() {
	// controls

	controls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	controls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	controls.addEventListener( 'change', cameraPositionChanged );

	controls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	controls.dampingFactor = 0.05;

	controls.enablePan = true;
	controls.enableZoom = true;

	controls.maxPolarAngle = Math.PI;
}

function cameraPositionChanged() {
	postStatus(`Camera position (${camera.position.x.toFixed(2)}, ${camera.position.y.toFixed(2)}, ${camera.position.z.toFixed(2)})`);
	
	// keep the raytracing sphere centred on the camera position
	// raytracingSphere.position.copy(camera.position.clone());	// TODO this doesn't seem to work as intended!?
}

async function toggleFullscreen() {
	if (!document.fullscreenElement) {
		document.documentElement.requestFullscreen().catch((err) => {
			postStatus(
				`Error attempting to enable fullscreen mode: ${err.message} (${err.name})`,
			);
		});
		// allow screen orientation changes
		// screen.orientation.unlock();
	} else {
		document.exitFullscreen();
	}
}

function showStoredPhoto() {
	gui.hide();
	renderer.domElement.style.visibility = "hidden";
	document.getElementById('takePhotoButton').style.visibility = "hidden";
	// document.getElementById('changePositionButton').style.visibility = "hidden";
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('backButton').style.visibility = "visible";
	document.getElementById('shareButton').style.visibility = "visible";
	document.getElementById('deleteButton').style.visibility = "visible";
	document.getElementById('storedPhoto').style.visibility = "visible";
	showingStoredPhoto = true;

	postStatus('Showing stored photo, '+storedPhotoDescription);
}

function showLivePhoto() {
	gui.show();
	renderer.domElement.style.visibility = "visible";
	document.getElementById('takePhotoButton').style.visibility = "visible";
	// document.getElementById('changePositionButton').style.visibility = "visible";
	if(storedPhoto) document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	document.getElementById('backButton').style.visibility = "hidden";
	document.getElementById('shareButton').style.visibility = "hidden";
	document.getElementById('deleteButton').style.visibility = "hidden";
	document.getElementById('storedPhoto').style.visibility = "hidden";
	showingStoredPhoto = false;

	postStatus('Showing live image');
}

function deleteStoredPhoto() {
	storedPhoto = null;

	showLivePhoto();

	postStatus('Stored photo deleted; showing live image');
}

function takePhoto() {
	try {
		click.play();

		storedPhoto = renderer.domElement.toDataURL('image/png');
		storedPhotoInfoString = getInfoString();

		storedPhotoDescription = name;
		// 
		document.getElementById('storedPhoto').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	
		postStatus('Photo taken; click thumbnail to view and share');
	} catch (error) {
		console.error('Error:', error);
	}	
}

async function share() {
	try {
		fetch(storedPhoto)
		.then(response => response.blob())
		.then(blob => {
			const file = new File([blob], name+storedPhotoDescription+'.png', { type: blob.type });

			// Use the Web Share API to share the screenshot
			if (navigator.share) {
				navigator.share({
					title: storedPhotoDescription,
					files: [file],
				});
			} else {
				postStatus('Sharing is not supported by this browser.');
			}	
		})
		.catch(error => {
			console.error('Error:', error);
			postStatus(`Error: ${error}`);
		});
	} catch (error) {
		console.error('Error:', error);
	}
}

/** 
 * Add a text field to the bottom left corner of the screen
 */
function createStatus() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	status.style.position = 'absolute';
	status.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	status.style.color = "White";
	status.style.fontFamily = "Arial";
	status.style.fontSize = "9pt";
	postStatus("Welcome!");
	status.style.bottom = 0 + 'px';
	status.style.left = 0 + 'px';
	status.style.zIndex = 1;
	document.body.appendChild(status);	
}

function postStatus(text) {
	status.innerHTML = '&nbsp;'+text;
	console.log('status: '+text);

	// show the text only for 3 seconds
	statusTime = new Date().getTime();
	setTimeout( () => { if(new Date().getTime() - statusTime > 2999) status.innerHTML = '&nbsp;'+name+', University of Glasgow, <a href="https://github.com/jkcuk/'+name+'">https://github.com/jkcuk/'+name+'</a>' }, 3000);
}

function getInfoString() {
	return `Lenslet array 1 (the closer array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;visible = ${raytracingSphereShaderMaterial.uniforms.visible1.value}<br>` +
		`&nbsp;&nbsp;period = ${raytracingSphereShaderMaterial.uniforms.period1.value.toFixed(4)}<br>` +
		`&nbsp;&nbsp;rotation angle = ${raytracingSphereShaderMaterial.uniforms.alpha1.value.toFixed(4)}&deg;<br>` +
		`&nbsp;&nbsp;focal length = ${raytracingSphereShaderMaterial.uniforms.focalLength1.value.toFixed(4)}<br>` +
		`&nbsp;&nbsp;radius = ${raytracingSphereShaderMaterial.uniforms.radius1.value.toFixed(4)}<br>` +
		`&nbsp;&nbsp;centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.x.toFixed(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.y.toFixed(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z.toFixed(4)})<br>` +
		`Lenslet array 2 (the farther array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;visible = ${raytracingSphereShaderMaterial.uniforms.visible2.value}<br>` +
		`&nbsp;&nbsp;period = ${raytracingSphereShaderMaterial.uniforms.period2.value.toFixed(4)} (&Delta;period = ${deltaPeriod.toFixed(4)})<br>` +
		`&nbsp;&nbsp;rotation angle = ${raytracingSphereShaderMaterial.uniforms.alpha2.value.toFixed(4)}&deg;<br>` +
		`&nbsp;&nbsp;focal length = ${raytracingSphereShaderMaterial.uniforms.focalLength2.value.toFixed(4)}<br>` +
		`&nbsp;&nbsp;radius = ${raytracingSphereShaderMaterial.uniforms.radius2.value.toFixed(4)}<br>` +
		`&nbsp;&nbsp;centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.x.toFixed(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.y.toFixed(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z.toFixed(4)}) (offset from confocal = ${offsetFromConfocal})<br>` +
		`camera position = (${camera.position.x.toFixed(4)}, ${camera.position.y.toFixed(4)}, ${camera.position.z.toFixed(4)})<br>` +
		`Horiz. FOV of user-facing(?) camera = ${fovVideoFeedU.toFixed(4)}&deg;<br>` +	// (user-facing) camera
		`Horiz. FOV of environment-facing(?) camera = ${fovVideoFeedE.toFixed(4)}&deg;<br>` +	// (environment-facing) camera
		`Horiz. FOV of screen = ${fovScreen.toFixed(4)}`;
}

function refreshInfo() {
	if(showingStoredPhoto) setInfo( storedPhotoInfoString );
	else setInfo( getInfoString() );

	if(info.style.visibility == "visible") setTimeout( refreshInfo , 100);	// refresh again a while
}

/** 
 * Add a text field to the top left corner of the screen
 */
function createInfo() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	info.style.position = 'absolute';
	info.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	info.style.color = "White";
	info.style.fontFamily = "Arial";
	info.style.fontSize = "9pt";
	setInfo("");
	info.style.top = 70 + 'px';
	info.style.left = 0 + 'px';
	info.style.zIndex = 1;
	document.body.appendChild(info);
	info.style.visibility = "hidden";
}

function setInfo(text) {
	info.innerHTML = text;
	console.log('info: '+text);
	// // show the text only for 3 seconds
	// infoTime = new Date().getTime();
	// setTimeout( () => { if(new Date().getTime() - infoTime > 2999) info.innerHTML = `` }, 3000);
	// info.style.visibility = "visible";
}

function toggleInfoVisibility() {
	switch(info.style.visibility) {
		case "visible":
			info.style.visibility = "hidden";
			break;
		case "hidden":
		default:
			info.style.visibility = "visible";
			refreshInfo();
	}
}