/* ============================================================================
 *  RUSLE SOIL EROSION MODEL  —  Google Earth Engine
 *  Revised Universal Soil Loss Equation:  A = R x K x LS x C x P
 * ----------------------------------------------------------------------------
 *  Author      : Furaha Albert Patrick  (Wanda Sphere Research)
 *  Region       : Rwanda (works for any FAO GAUL level-2 district worldwide)
 *  Resolution   : 10 m output (Sentinel-2 / ESA WorldCover scale)
 *  Licence      : MIT  —  free to use, adapt and share with attribution
 * ----------------------------------------------------------------------------
 *  HOW TO USE
 *    Change ONE line — the DISTRICT name below — and press Run.
 *    Optionally change COUNTRY and YEAR. Everything else is automatic.
 * ----------------------------------------------------------------------------
 *  SCIENTIFIC BASIS (why this version is not over-estimating)
 *    R  : Roose-type linear erosivity, R = 38.5 + 0.35 * P_annual
 *    K  : EPIC equation (Williams, 1995) from sand/silt/clay/organic carbon
 *    LS : McCool et al. (1987) SINE-based S factor + RUSLE beta slope-length.
 *         Sine (not tangent) prevents the classic over-estimation on the
 *         steep slopes that dominate Rwanda's "land of a thousand hills".
 *    C  : Durigon et al. (2014) tropical form, C = (1 - NDVI) / 2,
 *         computed on a FULL-YEAR Sentinel-2 composite (representative cover).
 *    P  : ESA WorldCover 10 m land cover + slope-dependent terracing factor.
 * ============================================================================ */


/* ----------------------------------------------------------------------------
 *  1. USER PARAMETERS  ——  edit these only
 * -------------------------------------------------------------------------- */
var COUNTRY      = 'Rwanda';      // FAO GAUL country name
var DISTRICT     = 'Muhanga';     // <<< CHANGE THIS to any district name
var YEAR         = 2024;          // analysis year
var EXPORT_SCALE = 10;            // metres; 10 = native, 30 = faster/lighter
var EXPORT_CRS   = 'EPSG:32735';  // UTM 35S (covers Rwanda); change if needed


/* ----------------------------------------------------------------------------
 *  2. DATA SOURCES
 * -------------------------------------------------------------------------- */
var CHIRPS = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD');
var DEM    = ee.Image('USGS/SRTMGL1_003');
var S2     = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var WCOVER = ee.ImageCollection('ESA/WorldCover/v200');
var SAND   = ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02');
var CLAY   = ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02');
var ORGC   = ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02');


/* ----------------------------------------------------------------------------
 *  3. STUDY AREA
 * -------------------------------------------------------------------------- */
var gaul = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level2');
var aoi  = gaul.filter(ee.Filter.and(
            ee.Filter.eq('ADM0_NAME', COUNTRY),
            ee.Filter.eq('ADM2_NAME', DISTRICT)));
var geom = aoi.geometry();

print('Study area:', DISTRICT + ', ' + COUNTRY);
print('Features found (should be >= 1):', aoi.size());
Map.centerObject(geom, 11);
Map.addLayer(aoi, {color: 'ffffff', fillColor: '00000000'}, DISTRICT + ' boundary');

var d1 = YEAR + '-01-01';
var d2 = (YEAR + 1) + '-01-01';


/* ----------------------------------------------------------------------------
 *  4. R FACTOR  —  Rainfall erosivity (MJ mm ha-1 h-1 yr-1)
 *     Near-uniform over one district = physically correct.
 * -------------------------------------------------------------------------- */
var P_annual = CHIRPS.filterDate(d1, d2).select('precipitation').sum().clip(geom);
var R = P_annual.multiply(0.35).add(38.5).rename('R');

Map.addLayer(R, {min: 380, max: 560,
  palette: ['fff7bc','fee391','fec44f','fe9929','d95f0e']}, 'R factor', false);


/* ----------------------------------------------------------------------------
 *  5. K FACTOR  —  Soil erodibility, EPIC equation (Williams 1995) -> SI
 * -------------------------------------------------------------------------- */
var SAN = SAND.select('b10').clip(geom).toFloat();                 // sand %
var CLA = CLAY.select('b10').clip(geom).toFloat();                 // clay %
var SIL = ee.Image(100).subtract(SAN).subtract(CLA);              // silt %
var OC  = ORGC.select('b10').clip(geom).toFloat().divide(10);     // g/kg -> %  (verify scaling)
var SN1 = ee.Image(1).subtract(SAN.divide(100));

var K = ee.Image().expression('(f1 * f2 * f3 * f4) * 0.1317', {
  'f1': ee.Image(0.2).add(ee.Image(0.3).multiply(
          SAN.multiply(-0.0256).multiply(ee.Image(1).subtract(SIL.divide(100))).exp())),
  'f2': SIL.divide(CLA.add(SIL)).pow(0.3),
  'f3': ee.Image(1).subtract(OC.multiply(0.25).divide(
          OC.add(OC.multiply(-2.95).add(3.72).exp()))),
  'f4': ee.Image(1).subtract(SN1.multiply(0.7).divide(
          SN1.add(SN1.multiply(22.9).add(-5.51).exp())))
}).rename('K').clip(geom);
K = K.max(0).min(0.07);

Map.addLayer(K, {min: 0.01, max: 0.045,
  palette: ['fff7bc','fee391','fec44f','fe9929','d95f0e','993404']}, 'K factor', false);


/* ----------------------------------------------------------------------------
 *  6. LS FACTOR  —  McCool et al. (1987) SINE-based topographic factor
 *     S = 10.8 sin(theta) + 0.03     (slope <  9 %)
 *     S = 16.8 sin(theta) - 0.50     (slope >= 9 %)
 *     L = (lambda / 22.13)^m ,  m = beta / (1 + beta)
 * -------------------------------------------------------------------------- */
var slopeDeg = ee.Terrain.slope(DEM.clip(geom));
var slopeRad = slopeDeg.multiply(Math.PI / 180);
var sinT     = slopeRad.sin();
var slopePct = slopeRad.tan().multiply(100).rename('s');

var S = ee.Image(0)
  .where(slopePct.lt(9),  sinT.multiply(10.8).add(0.03))
  .where(slopePct.gte(9), sinT.multiply(16.8).subtract(0.50));

var beta = sinT.divide(0.0896).divide(sinT.pow(0.8).multiply(3).add(0.56));
var m    = beta.divide(beta.add(1));
var L    = ee.Image(30).divide(22.13).pow(m);          // lambda = 30 m (SRTM cell)

var LS = L.multiply(S).rename('LS').clip(geom);
LS = LS.max(0).min(30);                                // safety cap for DEM spikes

Map.addLayer(LS, {min: 0, max: 12,
  palette: ['fff7bc','fee391','fec44f','fe9929','d95f0e','993404']}, 'LS factor', false);


/* ----------------------------------------------------------------------------
 *  7. C FACTOR  —  Durigon et al. (2014), full-year S2 composite
 * -------------------------------------------------------------------------- */
function maskS2(img) {
  var scl  = img.select('SCL');
  var good = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(7));
  return img.updateMask(good).divide(10000)
            .copyProperties(img, ['system:time_start']);
}
var s2col = S2.filterDate(d1, d2).filterBounds(geom)
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40)).map(maskS2);
print('Sentinel-2 scenes used:', s2col.size());

var ndvi = s2col.select(['B8','B4']).median()
                .normalizedDifference(['B8','B4']).rename('NDVI').clip(geom);
var C = ee.Image(1).subtract(ndvi).divide(2).rename('C').clip(geom).max(0).min(1);

Map.addLayer(ndvi, {min: 0, max: 0.85,
  palette: ['d73027','fdae61','fee08b','d9ef8b','1a9850']}, 'NDVI', false);
Map.addLayer(C, {min: 0, max: 0.7,
  palette: ['1a9850','a6d96a','fee08b','fdae61','d73027']}, 'C factor', false);


/* ----------------------------------------------------------------------------
 *  8. P FACTOR  —  ESA WorldCover 10 m + slope-dependent terracing
 * -------------------------------------------------------------------------- */
var lulc = WCOVER.mosaic().select('Map').rename('lulc').clip(geom);
var P = lulc.addBands(slopePct).expression(
  "(b('lulc') == 10) ? 0.80" +                            // tree cover
  ": (b('lulc') == 20) ? 0.85" +                          // shrubland
  ": (b('lulc') == 30) ? 0.90" +                          // grassland
  ": (b('s') <  2)  && (b('lulc') == 40) ? 0.55" +        // cropland by slope
  ": (b('s') <  5)  && (b('lulc') == 40) ? 0.50" +
  ": (b('s') <  8)  && (b('lulc') == 40) ? 0.50" +
  ": (b('s') < 12)  && (b('lulc') == 40) ? 0.60" +
  ": (b('s') < 16)  && (b('lulc') == 40) ? 0.70" +
  ": (b('s') < 20)  && (b('lulc') == 40) ? 0.80" +
  ": (b('s') >= 20) && (b('lulc') == 40) ? 0.90" +
  ": (b('lulc') == 50) ? 1.00" +                          // built-up
  ": (b('lulc') == 60) ? 1.00" +                          // bare / sparse
  ": (b('lulc') == 80) ? 0.00" +                          // water
  ": (b('lulc') == 90) ? 0.00" +                          // wetland
  ": 1.00"
).rename('P').clip(geom);

Map.addLayer(lulc, {min: 10, max: 100,
  palette: ['006400','ffbb22','ffff4c','f096ff','fa0000','b4b4b4',
            'f0f0f0','0064c8','0096a0','00cf75','fae6a0']}, 'WorldCover 10 m', false);
Map.addLayer(P, {min: 0, max: 1,
  palette: ['1a9850','a6d96a','fee08b','fdae61','d73027']}, 'P factor', false);


/* ----------------------------------------------------------------------------
 *  9. SOIL LOSS  —  A = R x K x LS x C x P
 * -------------------------------------------------------------------------- */
var A = R.multiply(K).multiply(LS).multiply(C).multiply(P).rename('Soil_Loss');

// Diagnostic — read before trusting the map. Healthy: median single digits-teens.
print('>>> Soil loss percentiles (t/ha/yr):', A.reduceRegion({
  reducer: ee.Reducer.percentile([10, 25, 50, 75, 90, 95, 99]),
  geometry: geom, scale: 90, maxPixels: 1e9, bestEffort: true}));
print('Mean soil loss (t/ha/yr):', A.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: geom,
  scale: 90, maxPixels: 1e9, bestEffort: true}).get('Soil_Loss'));

var PAL = ['2b83ba','91cba8','d9ef8b','ffffbf','fdae61','d7191c'];
Map.addLayer(A, {min: 0, max: 40, palette: PAL}, 'Soil loss (t/ha/yr)', false);


/* ----------------------------------------------------------------------------
 *  10. SEVERITY CLASSIFICATION  —  6 classes
 *      First break (10) = soil-loss tolerance (unsustainable above this).
 * -------------------------------------------------------------------------- */
var SLc = A.expression(
  "(b('Soil_Loss') <  5)  ? 1" +   // very slight
  ": (b('Soil_Loss') < 10) ? 2" +  // tolerable
  ": (b('Soil_Loss') < 20) ? 3" +  // moderate
  ": (b('Soil_Loss') < 40) ? 4" +  // high
  ": (b('Soil_Loss') < 80) ? 5" +  // severe
  ": 6"                            // extreme
).rename('SL_class').clip(geom);
Map.addLayer(SLc, {min: 1, max: 6, palette: PAL}, 'Soil loss class');


/* ----------------------------------------------------------------------------
 *  11. AREA PER CLASS  +  CHART
 * -------------------------------------------------------------------------- */
var areas = ee.Image.pixelArea().addBands(SLc).reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
  geometry: geom, scale: 90, maxPixels: 1e10, bestEffort: true});
var areaKm = ee.List(areas.get('groups')).map(function (it) {
  return ee.List([ee.Number(ee.Dictionary(it).get('sum')).divide(1e6).round()]); });
var classNames = ee.List(['Very slight (<5)','Tolerable (5-10)','Moderate (10-20)',
                          'High (20-40)','Severe (40-80)','Extreme (>80)']);

print(ui.Chart.array.values(areaKm, 0, classNames).setChartType('ColumnChart')
  .setOptions({title: 'Area per soil-loss class — ' + DISTRICT + ' (km2)',
    hAxis: {title: 'Severity class'}, vAxis: {title: 'Area (km2)'},
    legend: {position: 'none'}, colors: ['d7191c']}));


/* ----------------------------------------------------------------------------
 *  12. EXPORTS  (run from the Tasks tab)
 * -------------------------------------------------------------------------- */
Export.image.toDrive({
  image: A.toFloat(), description: DISTRICT + '_SoilLoss_' + YEAR,
  folder: 'earthengine', fileNamePrefix: DISTRICT + '_SoilLoss_' + YEAR,
  region: geom, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13});

Export.image.toDrive({
  image: SLc.toByte(), description: DISTRICT + '_SoilLossClass_' + YEAR,
  folder: 'earthengine', fileNamePrefix: DISTRICT + '_SoilLossClass_' + YEAR,
  region: geom, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13});


/* ----------------------------------------------------------------------------
 *  13. LEGEND
 * -------------------------------------------------------------------------- */
var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 15px'}});
legend.add(ui.Label({value: 'Soil loss (t ha-1 yr-1) — ' + DISTRICT + ' ' + YEAR,
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}}));
var row = function (c, n) {
  return ui.Panel({widgets: [
    ui.Label({style: {backgroundColor: '#' + c, padding: '8px', margin: '0 0 4px 0'}}),
    ui.Label({value: n, style: {margin: '0 0 4px 6px'}})],
    layout: ui.Panel.Layout.Flow('horizontal')}); };
var names = ['Very slight (<5)','Tolerable (5-10)','Moderate (10-20)',
             'High (20-40)','Severe (40-80)','Extreme (>80)'];
for (var i = 0; i < 6; i++) { legend.add(row(PAL[i], names[i])); }
Map.add(legend);

/* ============================================================================
 *  REFERENCES
 *  - Renard et al. (1997) RUSLE, USDA Agriculture Handbook 703.
 *  - McCool, D.K. et al. (1987) Revised slope steepness factor for the USLE.
 *      Transactions of the ASAE 30(5): 1387-1396.
 *  - Williams, J.R. (1995) The EPIC model (soil erodibility K).
 *  - Durigon, V.L. et al. (2014) NDVI-based C factor for tropical regions.
 *      International Journal of Remote Sensing 35(2): 441-453.
 *  - Karamage, F. et al. (2016) Extent of Cropland and Related Soil Erosion
 *      Risk in Rwanda. Int. J. Environ. Res. Public Health 13(5): 503.
 * ============================================================================ */
