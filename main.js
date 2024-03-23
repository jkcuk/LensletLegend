
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

let scene;
let aspectRatioVideoFeed = 4.0/3.0;
let renderer;
let videoFeed;
let camera;
let controls;
let raytracingSphereShaderMaterial;
	
// Nokia HR20, according to https://www.camerafv5.com/devices/manufacturers/hmd_global/nokia_xr20_ttg_0/
let fovVideoFeed = 68;	// (environment-facing) camera
let fovScreen = 68;

let cameraLensDistance = 10.0;
let raytracingSphereRadius = 20.0;
let offsetFromConfocal = 0.0;
let deltaPeriod = 0.0;

// the info text area
let info = document.createElement('div');
let infotime;	// the time the last info was posted

let gui;

// true if stored photo is showing
let showingStoredPhoto = false;

let storedPhoto;
let storedPhotoDescription;

// from https://github.com/4nt0nio/jpegcam
const click = new Audio('./shutter.mp3');

init();
animate();

function init() {
	// create the info element first so that any problems can be communicated
	createInfo();

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

	createVideoFeed();

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
	screen.orientation.addEventListener("change", recreateVideoFeed);

	addOrbitControls();

	// the controls menu
	createGUI();
}

function animate() {
	requestAnimationFrame( animate );

	if(!showingStoredPhoto) {
		// update uniforms
		updateUniforms();

		renderer.render( scene,  camera );
	}
}

function updateUniforms() {
	let tanHalfFovH, tanHalfFovV;
	if(aspectRatioVideoFeed > 1.0) {
		// horizontal orientation
		tanHalfFovH = Math.tan(0.5*fovVideoFeed*Math.PI/180.0);
		tanHalfFovV = Math.tan(0.5*fovVideoFeed*Math.PI/180.0)/aspectRatioVideoFeed;
	} else {
		// vertical orientation
		tanHalfFovH = Math.tan(0.5*fovVideoFeed*Math.PI/180.0)*aspectRatioVideoFeed;
		tanHalfFovV = Math.tan(0.5*fovVideoFeed*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.tanHalfFovH.value = tanHalfFovH;
	raytracingSphereShaderMaterial.uniforms.tanHalfFovV.value = tanHalfFovV;

	// calculate the separation between the two arrays, s = f1 + f2 + offsetFromConfocal
	let s = 
		raytracingSphereShaderMaterial.uniforms.focalLength1.value + 
		raytracingSphereShaderMaterial.uniforms.focalLength2.value +
		offsetFromConfocal;
	// arrange them symmetrically around z=0
	raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z = +0.5*s;
	raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z = -0.5*s - 0.001;

	// set the array periods
	raytracingSphereShaderMaterial.uniforms.period2.value = raytracingSphereShaderMaterial.uniforms.period1.value + deltaPeriod;
}

function createVideoFeed() {
	videoFeed = document.getElementById( 'videoFeed' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// environment-facing camera
		const constraints = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'environment'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraints ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeed.srcObject = stream;
			videoFeed.play();

			videoFeed.addEventListener("playing", () => {
				aspectRatioVideoFeed = videoFeed.videoWidth / videoFeed.videoHeight;
				updateUniforms();
				setInfo(`Camera resolution ${videoFeed.videoWidth} &times; ${videoFeed.videoHeight}`);
			});
		} ).catch( function ( error ) {
			setInfo(`Unable to access camera/webcam (Error: ${error})`);
		} );
	} else {
		setInfo( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}
}

/** create raytracing phere */
function addRaytracingSphere() {
	const videoFeedTexture = new THREE.VideoTexture( videoFeed );
	videoFeedTexture.colorSpace = THREE.SRGBColorSpace;

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
			focalLength2: { value: -1.0 },
			centreOfArray2: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			radius2: { value: 5.0 },	// radius of array 2
			videoFeedTexture: { value: videoFeedTexture }, 
			// cameraLensDistance: { value: cameraLensDistance },
			tanHalfFovH: { value: 1.0 },
			tanHalfFovV: { value: 1.0 }
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
			// uniform float cameraLensDistance;
			// uniform vec3 cameraPosition;
			
			// lenslet array 1
			uniform bool visible1;
			uniform float alpha1;	// rotation angle of array 1
			uniform float period1;	// period of array 1
			uniform float focalLength1;	// focal length of array 1
			uniform vec3 centreOfArray1;	// centre of array 1,  and principal point of lenslet (0, 0)
			uniform float radius1;	// radius of array 1

			// uniform float deltaZ12;	// z separation between arrays 1 and 2
			// uniform float offsetFromConfocal;	// offset from confocal configuration

			// lenslet array 2
			uniform bool visible2;
			uniform float alpha2;	// rotation angle of array 2
			uniform float period2;	// period of array 2
			uniform float focalLength2;	// focal length of array 2
			uniform vec3 centreOfArray2;	// centre of array 2,  and principal point of lenslet (0, 0)
			uniform float radius2;	// radius of array 2

			uniform sampler2D videoFeedTexture;
			uniform float tanHalfFovH;
			uniform float tanHalfFovV;
			
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
				// "normalised" version of d, scaled such that the z component is 1
				vec3 d1 = d/d.z;

				// calculate the intersection point with that array
				float deltaZ = centreOfArray.z - p.z;
				p = p + d1*deltaZ;

				// does the intersection point lie with in the radius?
				vec2 r = p.xy - centreOfArray.xy;
				float r2 = dot(r, r);	// length squared of vector r
				if(r2 < radius*radius) {
					// the intersection point lies inside the radius, so the component does something to the ray

					// which direction is the ray passing through it?
					if(d.z < 0.0) {
						// the ray is passing through array 1 in the direction for which it is designed; deal with it accordingly
	
						// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
						d = vec3(lensletArrayDeflect(-d1.xy, p.xy, centreOfArray.xy, alpha, period, focalLength), sign(d.z));
						// d = d.z*d1;

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} else {
						// the ray is passing through array 1 in the opposite direction for which it is designed
						b *= vec4(0.7, 0.3, 0.3, 1);
					}
				
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

				// is the first array visible?
				if(visible1) passThroughLensletArray(p, d, b, centreOfArray1,  radius1, alpha1, period1, focalLength1);

				// is the second array visible?
				if(visible2) passThroughLensletArray(p, d, b, centreOfArray2,  radius2, alpha2, period2, focalLength2);

				// does the ray intersect the (infinitely distant) camera image whose angular width and height is
				// given by arctan(2*tanHalfFovH) and arctan(2*tanHalfFovV?
				vec3 d1 = d/d.z;
				if((abs(d1.x) < tanHalfFovH) && (abs(d1.y) < tanHalfFovV)) {
					// yes, the ray intersects the image; take the pixel colour from the camera's video feed
					gl_FragColor = texture2D(videoFeedTexture, vec2(0.5-0.5*d1.x/tanHalfFovH, 0.5-0.5*d1.y/tanHalfFovV));
				} else {
					// no it doesn't; give the pixel a default colour, depending on which hemisphere it points towards
					if(d.z < 0.0) {
						gl_FragColor = vec4(0.53, 0.81, 0.92, 1.0);	// sky blue
					} else {
						gl_FragColor = vec4(1, 1, 0, 1.0);	// yellow
					}
				}
				gl_FragColor *= b;
			}
		`
	});
	const raytracingSphere = new THREE.Mesh( geometry, raytracingSphereShaderMaterial ); 
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
		'Point forward (in -<b>z</b> direction)': pointForward,
		'Field of view of screen (&deg;)': fovScreen,
		'Field of view of camera (&deg;)': fovVideoFeed,
		'Restart video streams': function() { 
			recreateVideoFeed(); 
			setInfo("Restarting video stream");
		}
	}

	const folderArray1 = gui.addFolder( 'Near lenslet array (in plane <i>z</i><sub>1</sub> = 0)' );

	folderArray1.add( params1, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible1.value = v; } );
	folderArray1.add( params1, 'Focal length, <i>f</i><sub>1</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.focalLength1.value = f; } );
	folderArray1.add( params1, 'Period, <i>p</i><sub>1</sub>', 0.01, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period1.value = p; } );
	folderArray1.add( params1, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha1.value = alpha/180.0*Math.PI; } );

	const folderArray2 = gui.addFolder( 'Far lenslet array' );

	folderArray2.add( params1, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible2.value = v; } );
	folderArray2.add( params2, 'Focal length, <i>f</i><sub>2</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.focalLength2.value = f; } );
	folderArray2.add( params2, '&Delta;<sub>period</sub>, <i>p</i><sub>2</sub> - <i>p</i><sub>1</sub>', -0.1, 0.1).onChange( (p) => { deltaPeriod = p; } );
	folderArray2.add( params2, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha2.value = alpha/180.0*Math.PI; } );
	folderArray2.add( params2, 'Offset from confocal', -0.25, 0.25).onChange( (o) => { offsetFromConfocal = o; } );


	const folderSettings = gui.addFolder( 'Other controls' );
	// folderSettings.add( params, 'Toggle show circles');
	// folderSettings.add( params, 'Swap arrays' );
	folderSettings.add( params, 'Point forward (in -<b>z</b> direction)');
	folderSettings.add( params, 'Field of view of screen (&deg;)', 10, 170, 1).onChange( setScreenFOV );   
	folderSettings.add( params, 'Field of view of camera (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeed = fov; updateUniforms(); });   
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
	setInfo('Pointing camera forwards (in -<b>z</b> direction)');
}

function onWindowResize() {
	screenChanged();
	setInfo(`window size ${window.innerWidth} &times; ${window.innerHeight}`);	// debug
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeed() {
	// stop current video stream...
	videoFeed.srcObject.getTracks().forEach(function(track) { track.stop(); });

	// ... and re-create a new one, hopefully of the appropriate size
	createVideoFeed();
}

function addOrbitControls() {
	// controls

	controls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	controls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	controls.addEventListener( 'change', () => { setInfo(`Camera position (${camera.position.x.toFixed(2)}, ${camera.position.y.toFixed(2)}, ${camera.position.z.toFixed(2)})`)} );

	controls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	controls.dampingFactor = 0.05;

	controls.enablePan = false;
	controls.enableZoom = false;

	// allowing control of the distance can result in the view being no longer 
	// centred on the origin, so don't allow it
	// controls.minDistance = cameraOutsideDistance;
	// controls.maxDistance = cameraOutsideDistance;

	controls.maxPolarAngle = Math.PI;
}

async function toggleFullscreen() {
	if (!document.fullscreenElement) {
		document.documentElement.requestFullscreen().catch((err) => {
			setInfo(
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

	setInfo('Showing stored photo, '+storedPhotoDescription);
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

	setInfo('Showing live image');
}

function deleteStoredPhoto() {
	storedPhoto = null;

	showLivePhoto();

	setInfo('Stored photo deleted; showing live image');
}

function takePhoto() {
	try {
		click.play();

		storedPhoto = renderer.domElement.toDataURL('image/png');

		storedPhotoDescription = `LensletLegend`;
		// 
		document.getElementById('storedPhoto').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	
		setInfo('Photo taken; click thumbnail to view and share');
	} catch (error) {
		console.error('Error:', error);
	}	
}

async function share() {
	try {
		fetch(storedPhoto)
		.then(response => response.blob())
		.then(blob => {
			const file = new File([blob], 'LensletLegend '+storedPhotoDescription+'.png', { type: blob.type });

			// Use the Web Share API to share the screenshot
			if (navigator.share) {
				navigator.share({
					title: storedPhotoDescription,
					files: [file],
				});
			} else {
				setInfo('Sharing is not supported by this browser.');
			}	
		})
		.catch(error => {
			console.error('Error:', error);
			setInfo(`Error: ${error}`);
		});
	} catch (error) {
		console.error('Error:', error);
	}
}

/** 
 * Add a text field to the bottom left corner of the screen
 */
function createInfo() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	info.style.position = 'absolute';
	info.style.backgroundColor = "rgba(0, 0, 0, 0.5)";	// semi-transparent white
	info.style.color = "White";
	info.style.fontFamily = "Arial";
	info.style.fontSize = "9pt";
	setInfo("Welcome!");
	info.style.bottom = 0 + 'px';
	info.style.left = 0 + 'px';
	info.style.zIndex = 1;
	document.body.appendChild(info);	
}

function setInfo(text) {
	info.innerHTML = text;
	console.log(text);

	// show the text only for 3 seconds
	infotime = new Date().getTime();
	setTimeout( () => { if(new Date().getTime() - infotime > 2999) info.innerHTML = `LensletLegend, University of Glasgow, <a href="https://github.com/jkcuk/LensletLegend">https://github.com/jkcuk/LensletLegend</a>` }, 3000);
}