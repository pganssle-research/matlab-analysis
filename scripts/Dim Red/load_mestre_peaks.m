function out = load_mestre_peaks(filename)
% Load all the mestrelab peaks output in XML file with the document style
% <root>
%		<DataSet file='' name='' description='' nDims=''>
%			<Peak>
%				<f dim='n'>peak</f>
%				<Value></Value>
%			</Peak>
%		</DataSet>
% </root>

if ~exist('filename', 'var') || ~exist(filename, 'file')
	filename = get_path(['mestre_hist.mat', {'*.xml'}]);
end

if ~exist(filename, 'file')
	error('Invalid file name.')
end

root = xmlread(filename);
removeIndentNodes(root.getChildNodes);
allDataSets = root.getElementsByTagName('DataSet');
nsets = allDataSets.getLength();

datasets = cell(nsets, 1);
i = 1;
for ii = 0:(nsets-1)
	ds = allDataSets.item(ii);
	
	attrs = ds.getAttributes();
	name = attrs.getNamedItem('name');
	if ~isempty(name)
		name = char(name.getValue());
	else
		name = '';
	end
	
	dfile = attrs.getNamedItem('file');
	if ~isempty(dfile)
		dfile = char(dfile.getValue());
	else
		dfile = '';
	end
	
	desc = attrs.getNamedItem('description');
	if ~isempty(desc)
		desc = char(desc.getValue());
	else
		desc = '';
	end
	
	nDims = attrs.getNamedItem('nDims'); % Only necessary one
	if isempty(nDims)
		nsets = nsets -1;
		continue;
	else
		nDims = str2double(nDims.getValue());
	end
	
	pks = ds.getChildNodes();
	nPeaks = pks.getLength();
	if nPeaks < 1
		nsets = nsets-1;
		continue;
	end
	
	peaks = zeros(nPeaks, nDims+1);
	dpos = ones(1, nDims+1)*NaN;
	j = 1;
	for jj = 0:(nPeaks-1)
		peak = pks.item(jj);
		
		if ~peak.hasChildNodes()
			nPeaks = nPeaks-1;
			continue;
		end
		
		cnodes = peak.getChildNodes();
		pos = dpos;
		for kk = 0:(cnodes.getLength()-1)
			node = cnodes.item(kk);
			
			if strcmp(node.getNodeName(), 'f')
				attrs = node.getAttributes();
				dim = attrs.getNamedItem('dim');
				if ~isempty(dim)
					dim = str2double(dim.getValue());
				else
					continue;
				end
				
				pos(dim) = str2double(node.getTextContent());
			elseif strcmp(node.getNodeName(), 'Value')
				pos(nDims+1) = str2double(node.getTextContent());
			end
		end
		
		if any(isnan(pos))
			nPeaks = nPeaks -1;
			continue
		end
		
		peaks(j, :) = pos;
		j = j+1;
	end
	
	if nPeaks == 0
		nsets = nsets-1;
		continue;
	elseif(size(peaks, 1) > nPeaks)
		peaks(nPeaks:end, :) = [];
	end
	
	data = struct('Name', name, 'File', dfile, 'Description', desc, ...
						'nDims', nDims, 'nPeaks', nPeaks, 'peaks', peaks);
	datasets{i} = data;
	i = i+1;
end
datasets(i:end) = [];

out = datasets;

function removeIndentNodes( childNodes )
	numNodes = childNodes.getLength;
	remList = [];
	for i = numNodes:-1:1
		theChild = childNodes.item(i-1);
		if (theChild.hasChildNodes)
			removeIndentNodes(theChild.getChildNodes);
		else
			if ( theChild.getNodeType == theChild.TEXT_NODE && ...
				  ~isempty(char(theChild.getData()))         && ...
				  all(isspace(char(theChild.getData()))))
				remList(end+1) = i-1; % java indexing
			end
		end
	end
	for i = 1:length(remList)
		childNodes.removeChild(childNodes.item(remList(i)));
	end

	
	
	