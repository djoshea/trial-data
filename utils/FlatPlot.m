classdef FlatPlot < matlab.graphics.chartcontainer.ChartContainer

    properties
        XProject matlab.lang.OnOffSwitchState = 'on'
        YProject matlab.lang.OnOffSwitchState = 'on'
        ZProject matlab.lang.OnOffSwitchState = 'on'

        XProjectLocation (1, 1) = "far" % "near", "far", or specific numeric location
        YProjectLocation (1, 1) = "far"
        ZProjectLocation (1, 1) = "far"

        ProjectChromaShift (1, 1) = 0; % max is about 0.2;
        ProjectLuminanceShift (1, 1) = 0;
        ProjectLuminanceMin (1, 1) = 0;
        ProjectChromaMax (1, 1) = 0.3;

        ProjectAlpha (1, 1) = 0.5
    end

    properties(SetAccess = private,Transient,NonCopyable)
        axh matlab.graphics.axis.Axes       
        axhOverlay matlab.graphics.axis.Axes     

        matrixTransform (:, :) double

        h_tracked (:, 1) matlab.graphics.primitive.Data
        
        xproj (1, 1) double
        xproj_epsilon (1, 1) double % xproj +/- some small epsilon towards the camera
        xproj_fn (1, 1) function_handle = @(v) v;
        xproj_behind_fn (1, 1) function_handle = @(v, childorderfrac) v;

        yproj (1, 1) double
        yproj_epsilon (1, 1) double % yproj +/- some small epsilon towards the camera
        yproj_fn (1, 1) function_handle = @(v) v;
        yproj_behind_fn (1, 1) function_handle = @(v, childorderfrac) v;

        zproj (1, 1) double
        zproj_epsilon (1, 1) double % zproj +/- some small epsilon towards the camera
        zproj_fn (1, 1) function_handle = @(v) v;
        zproj_behind_fn (1, 1) function_handle = @(v, childorderfrac) v;

        % projected objects
        h_xy (:, 1) matlab.graphics.primitive.Data
        h_yz (:, 1) matlab.graphics.primitive.Data
        h_xz (:, 1) matlab.graphics.primitive.Data

        hListeners
        installedCallbacks (1, 1) logical = false;
        currentlyRepositioningAxes (1, 1) logical = false;

        lastXLim (1, 2) double
        lastYLim (1, 2) double
        lastZLim (1, 2) double
    end

    properties(Dependent)
        ax
        axOverlay

        % Provide properties to support setting & getting
        XLimits (1,2) double
        XLimitsMode {mustBeMember(XLimitsMode,{'auto','manual'})}
        YLimits (1,2) double
        YLimitsMode {mustBeMember(YLimitsMode,{'auto','manual'})}
        ZLimits (1,2) double
        ZLimitsMode {mustBeMember(ZLimitsMode,{'auto','manual'})}
    end


     methods
%         function obj = FlatPlot(axh, varargin)
%             % Call superclass constructor method
%             obj@matlab.graphics.chartcontainer.ChartContainer(varargin{:});
% 
%             obj.axhOrig = axh;
%         end
    end

    methods(Access = protected)
        function setup(obj)
            % Get the layout and set the grid size
%             tcl = obj.getLayout();

            obj.axh = axes("Parent", gcf);
            hold(obj.axh, 'on');

            obj.axhOverlay = axes("Parent", gcf, 'Position', [0 0 1 1], 'HitTest', 'off');
            axis(obj.axhOverlay, "off");
            uistack(obj.axhOverlay, 'top');
            obj.axhOverlay.Color = "none";
            obj.axhOverlay.XLim = [0 1];
            obj.axhOverlay.YLim = [0 1];
        end
        
        function update(obj)
            %disp("Updating chart!")
            if ~isvalid(obj.axh)
                return;
            end
            % cache the matrix transform
            obj.matrixTransform = FlatPlot.axesTransformToFigure(obj.axh);

            obj.updateProjectLocations();

            % for each handle in h_tracked, generate projected shadows on the walls
            h_all = cat(1, obj.h_yz, obj.h_xz, obj.h_xy);
            if ~isempty(h_all)
                mask = isvalid(h_all);
                delete(h_all(mask));
            end
            
            nH = numel(obj.h_tracked);

            % figure out how many will be flattened behind
            behind = false(nH, 1);
            for iH = 1:nH
                behind(iH) = ismember(obj.h_tracked(iH).Type, "surface");
            end
            nBehind = nnz(behind);
            frac_behind = zeros(nH, 1);
            frac_behind(behind) = linspace(0, 1, nBehind);

            [h_yz, h_xz, h_xy] = deal(gobjects(nH, 1)); %#ok<*PROP> 
            for iH = 1:nH
                h = obj.h_tracked(iH);
                [h_yz(iH), h_xz(iH), h_xy(iH)] = obj.buildShadows(h, frac_behind(iH));
            end

            obj.h_yz = h_yz;
            obj.h_xz = h_xz;
            obj.h_xy = h_xy;
        end

        function updateProjectLocations(obj)
            % compute the locations of the projected data points
            axh = obj.axh;
            xl = axh.XLim;
            xmid = mean(xl);

            yl = axh.YLim;
            ymid = mean(yl);

            zl = axh.ZLim;
            zmid = mean(zl);
            
            % build matrix of points at the limits of each axis on the midpoint of the other two
            X = [xl, repmat(xmid, 1, 4)]';
            Y = [ymid, ymid, yl, ymid, ymid]';
            Z = [repmat(zmid, 1, 4), zl]';
            [~, ~, depth] = obj.ds2figFromTransform(obj.axh, obj.matrixTransform, X, Y, Z);
%             bigeps = 1e6 * eps;
%             bigeps = 0.4;
            lerp = @(v, slo, shi, dlo, dhi) (v - slo) ./ (shi - slo) .* (dhi - dlo) + dlo;


            % a small positive unit along each axis
            eps_gain = 1e-4;
            xepsilon = (xl(2) - xl(1)) * eps_gain;
            yepsilon = (yl(2) - yl(1)) * eps_gain;
            zepsilon = (zl(2) - zl(1)) * eps_gain;

            if depth(2) > depth(1)
                % xl(2) is closer
                [xnear, xfar] = deal(xl(2), xl(1));
                xtowards = xepsilon;
            else
                % xl(1) is closer
                [xnear, xfar] = deal(xl(1), xl(2));
                xtowards = -xepsilon;
            end
            if depth(4) > depth(3)
                % yl(2) is closer
                [ynear, yfar] = deal(yl(2), yl(1));
                ytowards = yepsilon;
            else
                % yl(1) is closer
                [ynear, yfar] = deal(yl(1), yl(2));
                ytowards = -yepsilon;
            end
            if depth(6) > depth(5)
                % zl(2) is closer
                [znear, zfar] = deal(zl(2), zl(1));
                ztowards = zepsilon;
            else
                % zl(1) is closer
                [znear, zfar] = deal(zl(1), zl(2));
                ztowards = -zepsilon;
            end

            % the projection functions below map coordinates within the axis limits into coordinates within epsilon of
            % the projected-onto axis wall
            % we also create proj_behind_fn which is used to project coordinates into a fixed slot based on the object child order

            % we shift the projection origin into the center of the axis a bit to make room for the projection range
            if isnumeric(obj.XProjectLocation)
                obj.xproj = obj.XProjectLocation;
            elseif strcmp(obj.XProjectLocation, "near")
                obj.xproj = xnear - xtowards;
            elseif strcmp(obj.XProjectLocation, "far")
                obj.xproj = xfar + 2.1*xtowards;
            else
                warning("Invalid value for XProjectLocation");
                obj.xproj = xfar + 2.1*xtowards;
            end
            
            % and then compress the entire axis range into 1 xepsilon of that range
            % while preserving the same ordering along the axis
            obj.xproj_epsilon = obj.xproj + xepsilon;
            obj.xproj_fn = @(x) lerp(x, xl(1), xl(2), obj.xproj, obj.xproj_epsilon);
            
            % this places everything behind the projected stuff above, while maintaining the childorder
            xproj_behind = obj.xproj - 1.1*xtowards;
            obj.xproj_behind_fn = @(x, childorderfrac) ones(size(x))*xproj_behind - (1-childorderfrac) * xtowards;
            
            if isnumeric(obj.YProjectLocation)
                obj.yproj = obj.YProjectLocation;
            elseif strcmp(obj.YProjectLocation, "near")
                obj.yproj = ynear - ytowards;
            elseif strcmp(obj.YProjectLocation, "far")
                obj.yproj = yfar + 2.1*ytowards;
            else
                warning("Invalid value for YProjectLocation");
                obj.yproj = yfar + 2.1*ytowards;
            end
            obj.yproj_epsilon = obj.yproj + yepsilon;
            obj.yproj_fn = @(y) lerp(y, yl(1), yl(2), obj.yproj, obj.yproj_epsilon);

            yproj_behind = obj.yproj - 1.1*ytowards;
            obj.yproj_behind_fn = @(y, childorderfrac) ones(size(y))*yproj_behind - (1-childorderfrac) * ytowards;

            if isnumeric(obj.ZProjectLocation)
                obj.zproj = obj.ZProjectLocation;
            elseif strcmp(obj.ZProjectLocation, "near")
                obj.zproj = znear - ztowards;
            elseif strcmp(obj.ZProjectLocation, "far")
                obj.zproj = zfar + 2.1*ztowards;
            else
                warning("Invalid value for ZProjectLocation");
                obj.zproj = zfar + 2.1*ztowards;
            end
            obj.zproj_epsilon = obj.zproj + zepsilon;
            obj.zproj_fn  = @(z) lerp(z, zl(1), zl(2), obj.zproj, obj.zproj_epsilon);

            zproj_behind = obj.zproj - 1.1*ztowards;
            obj.zproj_behind_fn = @(z, childorderfrac) ones(size(z))*zproj_behind - (1-childorderfrac) * ztowards;

            obj.lastXLim = xl;
            obj.lastYLim = yl;
            obj.lastZLim = zl;
        end

        function [h_yz, h_xz, h_xy] = buildShadows(obj, h, childorderfrac)
            h_yz = obj.buildShadow(h, "x", childorderfrac);
            h_xz = obj.buildShadow(h, "y", childorderfrac);
            h_xy = obj.buildShadow(h, "z", childorderfrac);
        end

        function hproj = buildShadow(obj, h, flattenAxis, childorderfrac)
            % this is the workhorse that draws specific objects on the side walls
            arguments
                obj
                h (1, 1) matlab.graphics.primitive.Data
                flattenAxis (1, 1) string {mustBeMember(flattenAxis, ["x", "y", "z"])} 
                childorderfrac (1, 1) double
            end

            switch h.Type
                case 'surface'
                    % special case, draw convex hull of projection as a filled patch

                    % project points
                    switch flattenAxis
                        case "x"
                            pts = [h.YData(:), h.ZData(:)];
                            inds = convhull(pts);

                            x = obj.xproj_behind_fn(h.XData(inds), childorderfrac);
                            y = pts(inds, 1);
                            z = pts(inds, 2);
                            
                        case "y"
                            pts = [h.XData(:), h.ZData(:)];
                            inds = convhull(pts);

                            x = pts(inds, 1);
                            y = obj.yproj_behind_fn(h.YData(inds), childorderfrac);
                            z = pts(inds, 2);
                        
                        case "z"
                            pts = [h.XData(:), h.YData(:)];
                            inds = convhull(pts);

                            x = pts(inds, 1);
                            y = pts(inds, 2);
                            z = obj.zproj_behind_fn(h.ZData(inds), childorderfrac);
                    end

                    hproj = fill3(x, y, z, h.FaceColor, Parent=h.Parent, FaceLighting='none', BackFaceLighting='reverselit');

                    % copy these props directly, we'll adjust below
                    props = ["FaceAlpha", "EdgeColor", "EdgeAlpha", "AmbientStrength", "DiffuseStrength", "SpecularStrength"];
                    for prop = props
                        hproj.(prop) = h.(prop);
                    end

                otherwise

                    hproj = copyobj(h, obj.axh);
                   
                    % project points
                    switch flattenAxis
                        case "x"
                            hproj.XData = obj.xproj_fn(h.XData);
        
                        case "y"
                            hproj.YData = obj.yproj_fn(h.YData);
                        
                        case "z"
                            hproj.ZData = obj.zproj_fn(h.ZData) + 0.1;
                    end

                    % alter appearance
                    adjustColor = @(color) FlatPlot.adjustColor(color, shiftL=obj.ProjectLuminanceShift, ...
                        shiftC=obj.ProjectChromaShift, ...
                        minL=obj.ProjectLuminanceMin, ...
                        maxC=obj.ProjectChromaMax, ...
                        alpha=obj.ProjectAlpha);
                    
                    props = ["Color", "CData", "FaceColor", "EdgeColor", "MarkerFaceColor", "MarkerEdgeColor"];
                    for prop = props
                        if isprop(hproj, prop)
                            hproj.(prop) = adjustColor(hproj.(prop));
                        end
                    end
            end

            hproj.Clipping = false;

            if isstruct(h.UserData) 
                if isfield(h.UserData, "shadow_props")
                    vals = h.UserData.shadow_props;
                    props = string(fieldnames(vals))';
                    for prop = props
                        if isprop(hproj, prop)
                            hproj.(prop) = vals.(prop);
                        end
                    end
                end

                if isfield(h.UserData, "shadow_props_" + flattenAxis)
                    vals = h.UserData.("shadow_props_" + flattenAxis);
                    props = string(fieldnames(vals))';
                    for prop = props
                        if isprop(hproj, prop)
                            hproj.(prop) = vals.(prop);
                        end
                    end
                end
            end

        end

        function addChild(obj, h)
            % adds a graphics handle to trackedChildren so that we can provide shadows for them
            obj.h_tracked = cat(1, obj.h_tracked, h);
        end
    end

    methods % Callback handling
        function installCallbacks(obj)
            % these work faster than listening on xlim and ylim, but can
            % not update depending on how the axis limits are set
            if isa(obj.axh, 'matlab.graphics.axis.Axes') 
%                 set(zoom(obj.axh),'ActionPreCallback',@obj.prePanZoomCallback);
%                 set(zoom(obj.axh),'ActionPostCallback',@obj.postPanZoomCallback);
            end
            
%             if ~isempty(obj.hListeners)
%                 delete(obj.hListeners);
%                 obj.hListeners = [];
%             end
%             
%             % listeners need to be cached so that we can delete them before
%             % saving.
%             hl(1) = addlistener(obj.axh, {'XDir', 'YDir', 'ZDir'}, 'PostSet', @obj.axisCallback);
%             hl(2) = addlistener(obj.axh, {'XLim', 'YLim', 'ZLim'}, 'PostSet', @obj.axisIfLimsChangedCallback);
%             obj.hListeners = hl;
            
            obj.installedCallbacks = true;
        end
        
        function uninstallCallbacks(obj)
            % remove all callbacks except cla listener
            
            % these work faster than listening on xlim and ylim, but can
            % not update depending on how the axis limits are set
            set(zoom(obj.axh),'ActionPreCallback',[]);
            set(zoom(obj.axh),'ActionPostCallback',[]);

            % delete axis limit and direction property listeners
            if ~isempty(obj.hListeners)
                delete(obj.hListeners);
                obj.hListeners = [];
            end

            obj.installedCallbacks = false;
        end

%         function prePanZoomCallback(obj, varargin)
%             % first, due to weird issues with panning, make sure we have
%             % the right auto axis for this update
%             if numel(varargin) >= 2 && isstruct(varargin{2}) && isfield(varargin{2}, 'Axes')
%                  axh = varargin{2}.Axes; %#ok<*PROPLC> 
%                  if obj.axh ~= axh
%                      return;
%                  end
%             end
%             obj.currentlyRepositioningAxes = true;
%             delete(obj.hListeners);
%             obj.hListeners = [];
%         end
% 
%         function postPanZoomCallback(obj, varargin)
%             % first, due to weird issues with panning, make sure we have
%             % the right auto axis for this update
%             if numel(varargin) >= 2 && isstruct(varargin{2}) && isfield(varargin{2}, 'Axes')
%                  axh = varargin{2}.Axes;
%                  if obj.axh ~= axh
%                      return;
%                  end
%             end
%             obj.currentlyRepositioningAxes = false;
%             obj.axisCallback(varargin{:});
%             obj.installCallbacks();
%         end 
        
        function axisIfLimsChangedCallback(obj, varargin)
            % similar to axis callback, but skips update if the limits
            % haven't changed since the last update
            if obj.isMultipleCall(), return, end
            
            if obj.currentlyRepositioningAxes
                % suppress updates when panning / zooming
                return;
            end
            
            % here we get clever. when panning or zooming, LocSetLimits is
            % used to set XLim, then YLim, which leads to two updates. We
            % check whether we're being called via LocSetLimits and then
            % don't update if we're setting the XLim, only letting the YLim
            % update pass through. This cuts our update time in half
            if numel(varargin) >= 1 && isa(varargin{1}, 'matlab.graphics.internal.GraphicsMetaProperty') 
                if strcmp(varargin{1}.Name, 'XLim')
                    % setting XLim, skip if in LocSetLimits
                    st = dbstack();
                    if ismember('LocSetLimits', {st.name})
                        return;
                    end
                end
            end

            if obj.checkLimsChanged()
                obj.axisCallback();
            end
        end

        function tf = checkLimsChanged(obj)
            tf = ~isequal(get(obj.axh, 'XLim'), obj.lastXLim) || ...
                 ~isequal(get(obj.axh, 'YLim'), obj.lastYLim) || ...
                 ~isequal(get(obj.axh, 'ZLim'), obj.lastZLim);
        end
        
        function axisCallback(obj, varargin)
            if obj.isMultipleCall(), return, end
            
             if numel(varargin) >= 2 && isstruct(varargin{2}) && isfield(varargin{2}, 'Axes')
                 axh = varargin{2}.Axes;
                 if obj.axh ~= axh
                     return;
                 end
             end
             obj.update();
        end
    end

    methods % client methods for addig to the plot
        function doUpdate(obj)
            obj.update();
        end

        function [Xf, Yf, Zf] = transform_ds2fig(obj, X, Y, Z)
            [Xf, Yf, Zf] = obj.ds2figFromTransform(obj.axh, obj.matrixTransform, X, Y, Z);
        end

        function h = plot3(obj, varargin)
            h = plot3(varargin{:}, Parent=obj.axh);
            obj.addChild(h);
        end

        function h = scatter3(obj, varargin)
            h = scatter3(varargin{:}, Parent=obj.axh);
            obj.addChild(h)
        end

        function h = scattercol(obj, X, varargin)
            h = obj.scatter3(X(:, 1), X(:, 2), X(:, 3), varargin{:});
            obj.addChild(h)
        end

        function h = covariance_ellipse(obj, X, args, surf_args)
            arguments
                obj
                X (:, 3) double
                args.confidence (1, 1) = 0.5;
                args.factor = [];
                args.n = 20;
                surf_args.?matlab.graphics.primitive.Surface
            end

            mu = mean(X, 1, 'omitnan');
            C = cov(X, 0, 'omitrows');
            
            % adapted from error_ellipse by AJ Johnson [ mathworks.com/matlabcentral/fileexchange/4705-error_ellipse ]          
            [eigvec,eigval] = eig(C);
            % Compute quantile for the desired percentile
            if isempty(args.factor)
                k = sqrt(chi2inv(args.confidence, 3)); % r is the number of dimensions (degrees of freedom)
            else
                k = args.factor;
            end

            [X,Y,Z] = ellipsoid(0,0,0,1,1,1, args.n);
            XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
            
            X(:) = k*XYZ(:,1)+mu(1);
            Y(:) = k*XYZ(:,2)+mu(2);
            Z(:) = k*XYZ(:,3)+mu(3);
            
            surf_args_c = namedargs2cell(surf_args);
            h = surface(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3, surf_args_c{:}, Parent=obj.axh);
            obj.addChild(h);
        end
    end

    methods % Pass-thru methods
        function ax = get.ax(obj)
            ax = obj.axh;
        end

        function ax = get.axOverlay(obj)
            ax = obj.axhOverlay;
        end

        function axis(obj, varargin)
            axis(obj.axh, varargin{:});
        end

        function box(obj, varargin)
            box(obj.axh, varargin{:});
        end

        function grid(obj, varargin)
            grid(obj.axh, varargin{:});
        end

        function varargout = view(obj, varargin)
            [varargout{1:nargout}] = view(obj.ax, varargin{:});
        end

        function varargout = xlim(obj,varargin)
            [varargout{1:nargout}] = xlim(obj.ax,varargin{:});
        end
        
        function varargout = ylim(obj,varargin)
            [varargout{1:nargout}] = ylim(obj.ax,varargin{:});
        end

        function varargout = zlim(obj,varargin)
            [varargout{1:nargout}] = zlim(obj.ax,varargin{:});
        end

        function varargout = xlabel(obj,varargin)
            [varargout{1:nargout}] = xlabel(obj.ax,varargin{:});
        end
        function varargout = ylabel(obj,varargin)
            [varargout{1:nargout}] = ylabel(obj.ax,varargin{:});
        end
        function varargout = zlabel(obj,varargin)
            [varargout{1:nargout}] = zlabel(obj.ax,varargin{:});
        end

        % set and get methods for XLimits and XLimitsMode
        function set.XLimits(obj, xlm)
            obj.ax.XLim = xlm;
        end
        function xlm = get.XLimits(obj)
            xlm = obj.ax.XLim;
        end
        function set.XLimitsMode(obj,xlmmode)
            obj.ax.XLimMode = xlmmode;
        end
        function xlm = get.XLimitsMode(obj)
            xlm = obj.ax.XLimMode;
        end
        
        % set and get methods for YLimits and YLimitsMode
        function set.YLimits(obj,ylm)
            obj.ax.YLim = ylm;
        end
        function ylm = get.YLimits(obj)
            ylm = obj.ax.YLim;
        end
        function set.YLimitsMode(obj, ylmmode)
            obj.ax.YLimMode = ylmmode;
        end
        function ylm = get.YLimitsMode(obj)
            ylm = obj.ax.YLimMode;
        end

        function zlm = get.ZLimits(obj)
            zlm = obj.ax.ZLim;
        end
        function set.ZLimitsMode(obj,zlmmode)
            obj.ax.YLimMode = zlmmode;
        end
        function zlm = get.ZLimitsMode(obj)
            zlm = obj.ax.ZLimMode;
        end
    end
    
    methods(Static)
        function [posNorm, posPixels, posCm] = getTrueAxesPosition(axh, outer, args)
            % based on plotboxpos by Kelly Kearney https://github.com/kakearney/plotboxpos-pkg
            %PLOTBOXPOS Returns the position of the plotted axis region
            %
            % pos = plotboxpos(h)
            %
            % This function returns the position of the plotted region of an axis,
            % which may differ from the actual axis position, depending on the axis
            % limits, data aspect ratio, and plot box aspect ratio.  The position is
            % returned in the same units as the those used to define the axis itself.
            % This function can only be used for a 2D plot.  
            %
            % Input variables:
            %
            %   h:      axis handle of a 2D axis (if ommitted, current axis is used).
            %
            % Output variables:
            %
            %   pos:    four-element position vector, in same units as h
            % Copyright 2010 Kelly Kearney
            % Check input

            arguments
                axh (1, 1) matlab.graphics.axis.Axes       
                outer = false; % true returns OuterPosition; false returns Position
                args.normRelativeToFigure = true; % false returns normalized to container, typically a tiled layout object
            end

            % Get position of axis in pixels
            currunits = axh.Units;
            axh.Units = 'Pixels';
            if outer
                axisPos = axh.OuterPosition;
            else
                axisPos = axh.Position;
            end
            
            % Calculate box position based axis limits and aspect ratios
            darismanual  = strcmpi(get(axh, 'DataAspectRatioMode'),    'manual');
            pbarismanual = strcmpi(get(axh, 'PlotBoxAspectRatioMode'), 'manual');
            if ~darismanual && ~pbarismanual
                % simple case
                posPixels = axisPos;
                axh.Units = 'normalized';
                posNorm = axh.Position;
                axh.Units = 'centimeters';
                posCm = axh.Position;
                axh.Units = currunits;
                return;
                
            else
                xlim = get(axh, 'XLim');
                ylim = get(axh, 'YLim');

                % Deal with axis limits auto-set via Inf/-Inf use

                if any(isinf([xlim ylim]))
                    hc = get(axh, 'Children');
                    hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
                    xdata = get(hc, 'XData');
                    if iscell(xdata)
                        xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
                        xdata = cat(1, xdata{:});
                    end
                    ydata = get(hc, 'YData');
                    if iscell(ydata)
                        ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
                        ydata = cat(1, ydata{:});
                    end
                    isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                                ~isinf(ydata) & ~isnan(ydata);
                    xdata = xdata(isplotted);
                    ydata = ydata(isplotted);
                    if isempty(xdata)
                        xdata = [0 1];
                    end
                    if isempty(ydata)
                        ydata = [0 1];
                    end
                    if isinf(xlim(1))
                        xlim(1) = min(xdata);
                    end
                    if isinf(xlim(2))
                        xlim(2) = max(xdata);
                    end
                    if isinf(ylim(1))
                        ylim(1) = min(ydata);
                    end
                    if isinf(ylim(2))
                        ylim(2) = max(ydata);
                    end
                end
                dx = diff(xlim);
                dy = diff(ylim);
                dar = get(axh, 'DataAspectRatio');
                pbar = get(axh, 'PlotBoxAspectRatio');
                limDarRatio = (dx/dar(1))/(dy/dar(2));
                pbarRatio = pbar(1)/pbar(2);
                axisRatio = axisPos(3)/axisPos(4);
                if darismanual
                    if limDarRatio > axisRatio
                        pos(1) = axisPos(1);
                        pos(3) = axisPos(3);
                        pos(4) = axisPos(3)/limDarRatio;
                        pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
                    else
                        pos(2) = axisPos(2);
                        pos(4) = axisPos(4);
                        pos(3) = axisPos(4) * limDarRatio;
                        pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
                    end
                elseif pbarismanual
                    if pbarRatio > axisRatio
                        pos(1) = axisPos(1);
                        pos(3) = axisPos(3);
                        pos(4) = axisPos(3)/pbarRatio;
                        pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
                    else
                        pos(2) = axisPos(2);
                        pos(4) = axisPos(4);
                        pos(3) = axisPos(4) * pbarRatio;
                        pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
                    end
                end
            end
            % Convert plot box position to the units used by the axis
            hparent = get(axh, 'parent');
            hfig = ancestor(hparent, 'figure'); % in case in panel or similar
            currax = get(hfig, 'currentaxes');

            temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hfig);
            posPixels = temp.Position;
            temp.Units = 'Normalized';
            posNorm = temp.Position;
            temp.Units = 'centimeters';

            if isa(hparent, 'matlab.graphics.layout.TiledChartLayout') && args.normRelativeToFigure
                % posNorm is in units normalized to the tiled chart layout's position, not the figure, which is what we
                % want
                tiled_pos = hparent.Position;
                posNorm = [tiled_pos(1) + posNorm(1), tiled_pos(2) + posNorm(2), tiled_pos(3) * posNorm(3), tiled_pos(4) * posNorm(4)];
            end

            posCm = temp.Position;
            delete(temp);
            axh.Units = currunits;
            set(hfig, 'currentaxes', currax);
        end

        function matrixTransform = axesTransformToFigure(axh)
            % this is a lightly modified version of AxesTransformToFigure
            % https://www.mathworks.com/matlabcentral/fileexchange/43896-3d-data-space-coordinates-to-normalized-figure-coordinates-conversion
            % Copyright (c) 2013, MinLong Kwong <komelong@gmail.com>
            
            arguments
                axh (1, 1) matlab.graphics.axis.Axes       
            end
           
            %%%% obtain data needed
            % camera
            viewAngle = get(axh, 'CameraViewAngle');
            viewTarget = get(axh, 'CameraTarget');
            viewPosition = get(axh, 'CameraPosition');
            viewUp = get(axh, 'CameraUpVector');
            % axes direction
            axesDirection = strcmp(get(axh, {'XDir', 'YDir', 'ZDir'}), 'normal');
            % data scale
            dataZLim = get(axh, 'ZLim');
            dataRatio = get(axh, 'DataAspectRatio');
            if any(dataRatio == 0), return, end
            plotBoxRatio = get(axh, 'PlotBoxAspectRatio');
            if any(plotBoxRatio == 0), return, end
            % is perspective
            isPerspective = strcmp(get(axh, 'Projection'), 'perspective');
            
            [positionNormal, positionPixel] = FlatPlot.getTrueAxesPosition(axh);
            
            % stretch
            stretchMode = strcmp(get(axh, {'CameraViewAngleMode', ...
                'DataAspectRatioMode', 'PlotBoxAspectRatioMode'}), 'auto');
            stretchToFill = all(stretchMode);
            stretchToFit = ~stretchToFill && stretchMode(1);
            stretchNone = ~stretchToFill && ~stretchToFit;

            %%%% model view matrix
            % move data space center to viewTarget point
            matrixTranslate = eye(4);
            matrixTranslate(1:3, 4) = -viewTarget;
            % rescale data
            % note: matlab will rescale data space by dividing DataAspectRatio
            %       and normalize z direction to 1 to makeup the 'PlotBox'
            scaleFactor = dataRatio / dataRatio(3) * (dataZLim(2) - dataZLim(1));
            scaleDirection = axesDirection * 2 - 1;
            matrixRescale = diag([scaleDirection ./ scaleFactor, 1]);
            % rotate the 'PlotBox'
            vecticesZUp = matrixRescale * ...
                [matrixTranslate * [viewPosition, 1]', [viewUp, 1]'];
            zVector = vecticesZUp(1:3, 1);
            upVector = vecticesZUp(1:3, 2);
            viewDistance = sqrt(dot(zVector, zVector));
            zDirection = zVector / viewDistance;
            yVector = upVector - zDirection * dot(zDirection, upVector);
            yDirection = yVector / sqrt(dot(yVector, yVector));
            matrixRotate = blkdiag( ...
                [cross(yDirection, zDirection), yDirection, zDirection]', 1);

            %%%% projection matrix
            % note: matlab will project the rotated 'PlotBox' to an area of 
            %       [-0.5, 0.5; -0.5, 0.5]
            matrixProjection = eye(4);
            matrixProjection(4, 3) = -isPerspective / viewDistance;
            projectionArea = 2 * tan(viewAngle * pi / 360) * viewDistance;
            matrixProjection = diag([ones(1, 3), projectionArea]) * matrixProjection;

            %%%% stretch matrix
            % stretch the projective 'PlotBox' into the position retangle of the axes
            % note: stretch will first detect data region
            if stretchToFill || stretchToFit
                plotBox = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]' - .5;
                plotBox = diag(plotBoxRatio / plotBoxRatio(3)) * plotBox;
                edgeVertices = matrixProjection * matrixRotate * [plotBox; ones(1, 8)];
                edgeVertices(1, :) = edgeVertices(1, :) ./ edgeVertices(4, :);
                edgeVertices(2, :) = edgeVertices(2, :) ./ edgeVertices(4, :);
                edgeVertices = edgeVertices(1:2, :)';
                % note: the low boundary and the high boundary of data region may be
                %       difference in perspective projection, so the figure should move
                %       to center, but here no need to do so, because matlab ignore it
                dataRegion = max(edgeVertices) - min(edgeVertices);
                % note: matlab have a strange addition stretch in stretch to fit mode.
                %       one side of the data region will hit the position rectangle,
                %       and matlab will assume data region of that side to be 1 keeping
                %       aspect ratio.
                if stretchToFit
                    strangeFactor = dataRegion ./ positionPixel(3:4);
                    if strangeFactor(1) > strangeFactor(2)
                        dataRegion = dataRegion / dataRegion(1);
                    else
                        dataRegion = dataRegion / dataRegion(2);
                    end
                end
            else
                % note: if no stretch, it will use projection area as data region
                dataRegion = [1, 1];
            end
            % note: stretch than apply a stretchFactor to the data, such that it fit
            %       in the axes position retangle
            if stretchToFit || stretchNone
                stretchFactor = dataRegion ./ positionPixel(3:4);
                stretchFactor = stretchFactor / max(stretchFactor);
            else
                stretchFactor = [1, 1];
            end
            matrixStretch = diag([stretchFactor ./ dataRegion, 1, 1]);

            %%%% view port matrix
            matrixViewPort = diag([positionNormal(3:4), 1, 1]);
            matrixViewPort(1:2, 4) = positionNormal(1:2) + positionNormal(3:4) / 2;

            %%%% return transformation matrix
            matrixTransform = matrixViewPort * matrixStretch * matrixProjection * ...
                matrixRotate * matrixRescale * matrixTranslate;
        end
    
        function coordLiner = adjustLogScaleData(axh, coordLogOrLiner)
            % this is a lightly modified version of AdjustLogScaleData
            % https://www.mathworks.com/matlabcentral/fileexchange/43896-3d-data-space-coordinates-to-normalized-figure-coordinates-conversion
            % Copyright (c) 2013, MinLong Kwong <komelong@gmail.com>

            arguments
                axh (1, 1) matlab.graphics.axis.Axes       
                coordLogOrLiner (:, :, :) double
            end
            
            % adjust axes data for log scale
            coordLiner = coordLogOrLiner;
            % obtain data
            isLogScale = strcmp(get(axh, {'XScale', 'YScale', 'ZScale'}), 'log'); 
            lim = get(axh, {'XLim', 'YLim', 'ZLim'});
            % adjust
            for i = 1:3
                if isLogScale(i)
                    if lim{i}(1) <= 0
                        error(['low boundary of ''', 'W'+i, 'Lim'' should not less', ...
                            ' than or equal to 0 in log scale axes, consider set ''', ...
                            'W'+i, 'LimMode'' to ''auto''.']);
                    end
                    rate = (log(coordLogOrLiner(i, :)) - log(lim{i}(1))) ...
                        / (log(lim{i}(2)) - log(lim{i}(1)));
                    coordLiner(i, :) = rate * (lim{i}(2) - lim{i}(1)) + lim{i}(1);
                end
            end
        end

        function [Xf, Yf, Zf] = ds2figFromTransform(axh, matrixTransform, X, Y, Z)
            % dataCoords is nPts x 3
            % this is a lightly modified version of ds2fig
            % https://www.mathworks.com/matlabcentral/fileexchange/43896-3d-data-space-coordinates-to-normalized-figure-coordinates-conversion
            % Copyright (c) 2013, MinLong Kwong <komelong@gmail.com>
            
            arguments
                axh (1, 1) matlab.graphics.axis.Axes      
                matrixTransform (:, :) double
                X double
                Y double
                Z double
            end
            
            % get data space coordinates from arguments
            if isempty(Z)
                Z = zeros(size(X));
            end

            % check if all arguments are same size
            dataSize = size(X);
            if ~(all(size(Y) == dataSize) && all(size(Z) == dataSize))
                error('X,Y,Z should be the same size');
            end
            
            % build data space coordinates
            X = X(:);
            Y = Y(:);
            Z = Z(:);
            dataCoord = [X, Y, Z, ones(size(X))]';
            dataCoord = FlatPlot.adjustLogScaleData(axh, dataCoord);

            % transform and build figure coordinates
            figureCoord = matrixTransform * dataCoord;
            % perspective division
            Xf = reshape(figureCoord(1, :) ./ figureCoord(4, :), dataSize);
            Yf = reshape(figureCoord(2, :) ./ figureCoord(4, :), dataSize);
            Zf = reshape(figureCoord(3, :) ./ figureCoord(4, :), dataSize);
        end

        function rgb = adjustColor(rgb, lch_args, args)
            arguments
                rgb 

                lch_args.shiftL = 0;
                lch_args.gainL = 1;
                lch_args.minL = 0;
                lch_args.maxL = 1;
        
                lch_args.shiftC = 0;
                lch_args.gainC = 1;
                lch_args.minC = 0;
                lch_args.maxC = 0.35;
                
                lch_args.shiftH = 0;
                lch_args.minH = 0;
                lch_args.maxH = 360;

                args.alpha = 1;
            end
            if ~isnumeric(rgb) || isempty(rgb)
                return;
            end

            lch_args_c = namedargs2cell(lch_args);
            rgb = TrialDataUtilities.Color.adjust_oklch(rgb, lch_args_c{:});

            rgb = TrialDataUtilities.Color.computeWithAlpha(rgb, args.alpha);
        end

        function flag = isMultipleCall()
            % determine whether callback is being called within itself
            flag = false; 
            % Get the stack
            s = dbstack();
            if numel(s) <= 2
                % Stack too short for a multiple call
                return
            end

            % How many calls to the calling function are in the stack?
            names = {s(:).name};
            TF = strcmp(s(2).name,names);
            count = sum(TF);
            if count>1
                % More than 1
                flag = true; 
            end
        end
    end
end


