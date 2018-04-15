// Task --------------------------------------------------

const R_a = 20.0 * Math.pow(10.0, -2);  // 20cm - interior radius
const R_b = 100.0 * Math.pow(10.0, -2); // 100cm - external radius
const G = 1.5 * Math.pow(10.0, 5);      // 1.5*10^5 W/m^3
const T_b = 20.0;                       // 20C
const Q_a = 2500.0;                     // 2500 W/m^2
const K = 50.0;                         // W/(m*deg)

const p = 2.0;                          // formula's coefficient
const N = 50;                           // amount of parts
const delta = (R_b - R_a) / N;          // amount of iterations

//--------------------------------------------------------

// T(r)-? | Q(R_a)-? | Q(R_b)-?

// Solution ----------------------------------------------

var a = [];
var b = [];
var c = [];
var d = [];
var T = [];

createGeometry(a, b, c, d);

Thomas(N + 1, a, b, c, d, T);

outputStreamDOM(calculateStreamInterior(T, K, delta), "Q(R_a)");
outputStreamDOM(calculateStreamExternal(T, K, delta), "Q(R_b)");


// output in Google Charts TM
var resultX = new matrixArray(1, 0);
var resultY = new matrixArray(1, 0);

for (let i = 0; i < N + 1; i++) {
    resultX[0][i] = R_a + i * delta;
    resultY[0][i] = T[i];
}

//--------------------------------------------------------



/** Form a nodal point system of N equal parts -----------
  * Where diagonals (vectors) are stored as follows:
  * a = [0; 49] -- bottom,
  * b = [0; 50] -- middle,
  * c = [0; 49] -- top,
  * d = [0; 50] -- right side.
  */

function createGeometry(a, b, c, d) {
    var rhsValue = - G * Math.pow(delta, 2) / K;

    for (let i = 1; i < N; i++) {
        var coef = (p * delta / 2.0) / (R_a + delta * i);
        a[i - 1] = 1 - (p * delta / 2.0) / (R_a + delta * i);
        b[i] = - 2.0;
        c[i] = 1 + coef;
        d[i] = rhsValue;
    }

    // left border -- K * dT/dr | {r = Rb} = - Q_a
    b[0] = -2.0;
    c[0] = 2.0;
    var coef = (p * delta / 2.0) / (R_a + 0.0 * delta);
    d[0] = rhsValue - (2 * Q_a * delta / K) * (1 - coef);

    // right border -- T(Rb) = T_b;
    a[N - 1] = 0.0;
    b[N] = 1.0;
    d[N] = T_b;
}



// Thomas algorithm for 3-diagonal matrices --------------

function Thomas(N, a, b, c, d, T) {
    let m = 0.0;
    for (let i = 1; i < N; i++) {
        m = a[i - 1] / b[i - 1];
        b[i] = b[i] - m * c[i - 1];
        d[i] = d[i] - m * d[i - 1];
    }

    T[N - 1] = d[N - 1] / b[N - 1];

    for (let i = N - 2; i >= 0; i--) {
        T[i] = (d[i] - c[i] * T[i + 1]) / b[i];
    }
}



// Utility functions -------------------------------------

// temperature stream calculation using derivative's decomposition
function calculateStreamInterior(T, K, delta) {
    // 3 points: m = 2.0, j = 0; A0 = -3, A1 = 4, A2 = -1
    // return streamInterior = - K*(-3.0*T[0] + 4.0*T[1] - 1.0*T[2]) / (2.0*delta);

    // 6 points: ...
    return streamInterior = -K*(-274.0*T[0] + 600.0*T[1] -
        600.0*T[2] + 400.0*T[3] - 150.0*T[4] + 24.0*T[5]) / (120.0*delta);
}

function calculateStreamExternal() {
    // 3 points: m = 2.0, j = 2; A0 = 1, A1 = -4, A2 = 3
    // return streamExternal = K*(1.0*T[N-2] - 4.0*T[N-1] + 3.0*T[N]) / (2.0*delta);

    // 6 points: ...
    return streamExternal = K*(-24.0*T[N - 5] + 150.0*T[N - 4] -
        400.0*T[N - 3] + 600.0*T[N-2]-600.0*T[N-1]+274.0*T[N]) / (120.0*delta);
}

// matrix-like arrays used in google charts
function matrixArray(rows, columns) {
    var arr = new Array();
    for(var i = 0; i < rows; i++)
        arr[i] = new Array();
    return arr;
}

// object info's dump for debug (similar to PHP) 
function varDump(obj, string) {
    var out = '';
    for (var i in obj)
        out += "" + string + "\t" + i + ":\t\t " + obj[i] + "\n";

    console.log('Object: \n' + out);
}

// output of a stream into the HTML
function outputStreamDOM(stream, title) {
    var outputString = "Stream " + title + " : " + stream / Math.pow(10.0, 3) + " [kW / m^2]";

    console.log(outputString);

    var streamInteriorDOM = document.createElement('h2');
    streamInteriorDOM.innerHTML = outputString + '<br />';
    document.body.appendChild(streamInteriorDOM);
}



// Tables ------------------------------------------------

google.charts.load('current', {'packages':['table']});
google.charts.setOnLoadCallback(drawTable);

function drawTable() {

    data = new google.visualization.DataTable();

    data.addColumn('string', 'Cross-Section of the Sphere [m]');
    data.addColumn('string', 'Temperature Field [°С]');
    data.addColumn('number', 'Step');

    for(let i = 0; i < N + 1; i++)
        data.addRows([ ['' + resultX[0][i], '' + resultY[0][i], i] ]);

    var tempName = 'table_div'
    table = new google.visualization.Table(document.getElementById(tempName));

    table.draw(data, {showRowNumber: false, width: '100%', height: '100%'});
}



// Graphs ------------------------------------------------

google.charts.load('current', {packages: ['corechart', 'line']});
google.charts.setOnLoadCallback(drawBackgroundColor);

function drawBackgroundColor() {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'X');
    data.addColumn('number', '1');

    for(let i = 0; i < N + 1; i++)
        data.addRows([ [resultX[0][i],  resultY[0][i]] ]);

    var options = {
    hAxis: {
        title: 'r [ m ]',
        minValue: R_a,
        // maxValue: R_b,
        format: 'decimal',
        gridlines: { count: 5 },
    },
    vAxis: {
        title: 'T [ °С ]',
        minValue: 0,
        format: 'decimal', // scientific
        gridlines: { count: 10 }
    },
    backgroundColor: 'white',
    };

    var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));
    chart.draw(data, options);
}
