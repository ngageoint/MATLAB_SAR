# MATLAB SAR Toolbox
The MATLAB SAR Toolbox is a basic MATLAB library to read, write, display, and do simple processing of complex SAR data using the NGA [SICD](http://www.gwg.nga.mil/ntb/baseline/docs/SICD/) format.  It has been released by NGA to encourage the use of SAR data standards throughout the international SAR community.  The MATLAB SAR Toolbox complements the [SIX library](https://github.com/ngageoint/six-library) (C++) and [SarPy](https://github.com/ngageoint/sarpy) (Python), which are implemented in other languages but have similar goals.

Some sample SICD files are available [here](https://github.com/ngageoint/six-library/wiki/Sample-SICDs), although since the MATLAB SAR Toolbox can make its own SICDs or read other formats with an API that makes them look as if they were SICD, one doesn't really need an actual SICD file to use this toolbox and benefit from the SICD metadata standard.

A sampling of some of the functionality available in the toolbox is provided below.

### File I/O
The toolbox has readers and SICD converters for the following complex SAR formats:
* ALOS PALSAR-2
* COSMO-SkyMed
* Complex NITF
* GFF (Sandia)
* NISAR
* Sentinel-1
* SICD
* SIO
* RADARSAT-2
* RADARSAT Constellation Mission (RCM)
* TerraSAR-X

To read data from a file:
```
readerobj = open_reader(filename);  % Toolbox will detect file format and parse appropriately
sicdmetadata = readerobj.get_meta();  % Metadata from all formats will be converted to SICD metadata structure
complex_data = readerobj.read_chip([first_column last_column], [first_row last_row]);
readerobj.close();
% Or you can use read_complex_data.m in a single line of code if you won't be reading multiple chunks from an image.
```
It is important to stress that, since the metadata returned will always be in SICD format regardless of the original format of the data, as long as one writes code to the SICD standard, that code will generically process all of the above formats.

To write data to SICD:
```
writer_object = NITFWriter(filename, sicdmetadata);
writer_object.write_chip(data, [first_row first_column]);
% The file will be closed when writer_object is deleted or goes out of scope.
```

Note also that convert_complex_data.m will take a complex image from any format recognized and convert it to SICD file(s).

### SICD Validation
For those wanting to produce SICD files from their own SAR data sources, validate_sicd.m is a function that runs a set of over one hundred tests to check for validity in SICD files.  It does this through a set of checks for internal consistency with SICD metadata, as well as a set of interactive tests for testing consistency between the SICD metadata and the SICD pixel data.

### Processing
Geometry/Projections provides implementations of precise scene-to-image and image-to-scene projections using the SICD sensor model.

Examples using interactive tools are provided for subaperture processing (ApertureTool.m), radiometric measurement (RCSTool.m), and polarimetric visualization (PolTool.m).

With all of these tools, every effort was made to assure they were generic to a wide variety of SAR data types (spotlight/stripmap, zero Doppler range migration/polar format algorithm, etc.)

### Visualization
The Taser tool (or TaserClean, if you prefer) provides a way to easily and quickly browse/zoom/pan through large complex SAR datasets that might be too large and cumbersome to fit into memory and then potentially call other tools on selections of that data.  The component in those tools that enables this is a reusable component (hg_mitm_viewer.m) that can be easily added into new tools (as was done with the RCSTool and PolTool).

There are also tools for visualizing a SAR collection geometry in KML (Image2KMLGUI.m), browsing image metadata (MetaView.m), and generating a "metaicon" for summarizing a SAR collection (MIM.m).

### Origin
The MATLAB SAR Toolbox was developed at the National Geospatial-Intelligence Agency (NGA). The software use, modification, and distribution rights are stipulated within the MIT license.

### Pull Requests

If you'd like to contribute to this project, please make a pull request. We'll review the pull request and discuss the changes. All pull request contributions to this project will be released under the MIT license.

Software source code previously released under an open source license and then modified by NGA staff is considered a "joint work" (see 17 USC § 101); it is partially copyrighted, partially public domain, and as a whole is protected by the copyrights of the non-government authors and must be released according to the terms of the original open source license.
