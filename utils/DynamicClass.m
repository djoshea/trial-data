classdef (HandleCompatible) DynamicClass
% This class is very similar to dynamicprops, in that it allows classes or instances
% to respond to specific properties that are not defined in the class definition.
% The primary difference with dynamicprops is that it does not inherit handle,
% which allows it to be used on value classes by inheritance. 
% It handles a variety of different .property, .method
% () indexing, and {} indexing initiated from outside this class
% in a way that correctly

    properties(Constant, Hidden)
        % each of the methods below should return this to indicate that a particular
        % operation is not supported or a dynamic property/method name is not
        % recognized
        NotSupported = DynamicClassEnumeration.NotSupported;
    end
    
    properties(Hidden)
        dynamicClassActive = true;
    end

    % Override if you'd like dynamic access to take precedence over declared properties/methods
    methods(Hidden, Access=protected) 
        % allow mapDynamicPropertyAccess to handle accessing of declared properties
        % allow mapDynamicPropertyAssign to handle assiging declared properties
        % if false, declared properties will always be accessed/assigned directly
        function tf = allowDynamicPropertyOverride(obj) %#ok<MANU>
            tf = false;
        end

        % allow mapDynamicMethod to handle calling of declared methods
        % if false, declared methods will always be called directly
        function tf = allowDynamicMethodOverride(obj) %#ok<MANU>
            tf = false;
        end
    end

    % Override as many as you'd like for various forms of dynamic access
    methods(Hidden) % Methods for accessing (via subsref)
        % You should override any or all of these methods as needed
        % The implementations below simply return DynamicClass.NotSupported
       
        % given a dynamic property name and a set of idx, return the value of 
        % this property, or return DynamicClass.NotSupported if name is not
        % recognized as a valid property name
        %
        % indexingApplied indicates whether the property value
        % has already been indexed using idx, or whether the idx were ignored.
        %
        % With indexingApplied == false:
        %   val = obj.name(idx) --> val = obj.mapDynamicPropertyAccess(name); val = val(idx) or val = val{idx};
        % With includeIndexing == true:
        %   val = obj.name(idx) --> val = obj.mapDynamicPropertyAccess(name, idx);
        %
        % The purpose of this is that is often much faster to request only the 
        % requested idx of the dynamic property rather than generating all of them
        % and then selecting the idx to return. includeIndexing indicates whether
        % the method has opted to use this ability itself or ignore the idx and
        % defer indexing to our subsref implementation
        function [value, appliedNext] = mapDynamicPropertyAccess(obj, name, typeNext, subsNext) %#ok<INUSD>
            value = DynamicClass.NotSupported;
            appliedNext = false;
        end

        % return a function that should be called for this dynamic method name
        % or return DynamicClass.NotSupported if the method name is not recognized
        %
        % Note that the function will receive the name as the first argument,
        % followed by the rest of the arglist, i.e.
        %
        % varargout{:} = obj.name(varargin{:}) --> varargout{:} = fn(name, varargin{:})
        %
        function fn = mapDynamicMethod(obj, name) %#ok<INUSD>
            fn = DynamicClass.NotSupported;
        end

        % handle parenthetical indexing of the base class
        % obj(idx1, idx2, ...) --> obj.parenIndex({idx1, idx2, ...})
        % note that you should probably also implement size(), length(), and end() 
        % If parenthetical indexing is not supported, return DynamicClass.NotSupported
        %
        % typeNext and subsNext will typically contain either '', {} (in which case
        % you should ignore them and set appliedNext = false, or will have the values:
        % typeNext = '.' and subsNext = someFieldName. If you would like to process
        % the field reference (e.g. obj(idx).field) for efficiency, use the list of
        % subscripts {idx1, idx2, ...} in subsNext and set appliedNext = true;
        function [result, appliedNext] = parenIndex(obj, subs, typeNext, subsNext) %#ok<INUSD>
            result = DynamicClass.NotSupported;
            appliedNext = false;
        end

        % handle curly brace indexing of the base class
        % obj{idx1, idx2, ...} --> obj.cellIndex({idx1, idx2, ...})
        % If cell array indexing is not supported, return DynamicClass.NotSupported
        %
        % You must fill one entry of resultCell for each index requested (i.e. 
        % the size must be as specified by subs), otherwise subsref will throw
        % an error
        function [result, appliedNext] = cellIndex(obj, subs, typeNext, subsNext) %#ok<INUSD>
            result = DynamicClass.NotSupported;
            appliedNext = false;
        end
    end

    % Override as many as you'd like for various forms of dynamic assignment
    methods(Hidden) % Methods for assignment (via subsasgn)
        function obj = dynamicPropertyAssign(obj, name, value, s) %#ok<INUSD>
            obj = DynamicClass.NotSupported;
        end

        function obj = parenAssign(obj, subs, value, s) %#ok<INUSD>
            obj = DynamicClass.NotSupported;
        end

        function obj = cellAssign(obj, subs, value, s) %#ok<INUSD>
            obj = DynamicClass.NotSupported;
        end
    end

    % Internal Implementation
    methods(Sealed)
        function [result, s] = mapDeclaredPropertyAccess(obj, name, meta, s)
            propInfo = meta.PropertyList;
            idx = find(strcmp(name, {propInfo.Name}));
            assert(~isempty(idx), 'Could not find property info for %s', name);
            if strcmp(propInfo(idx).GetAccess, 'public')
                % public property, simply access it and store as intermediate
                result = obj.(name);
                s = s(2:end);
            else
                error('Property %s is not publicly accessible', name);
            end
        end

        function [result, s, returnImmediately] = mapDeclaredMethodCall(obj, name, meta, s, nargout)
            returnImmediately = false;
            methodInfo = meta.MethodList;
            idx = find(strcmp(name, {methodInfo.Name}));
            assert(~isempty(idx), 'Could not find property info for %s', name);

            if strcmp(methodInfo(idx).Access, 'public')
                % public method, call it and store as intermediate
                if length(s) > 1 && strcmp(s(2).type, '()')
                    % method call with argument list
                    args = s(2).subs;
                    s = s(3:end);
                else
                    % method call without arguments
                    args = {};
                    s = s(2:end);
                end
                
                % does further subsref indexing happen after this?
                if isempty(s)
                    % no more after this, call the method and return
                    % this syntax calls the method and preserves
                    % the correct nargout
                    [result{1:nargout}] = obj.(name)(args{:});
                    returnImmediately = true;
                    return;
                else
                    % more subsref indexing happens after this, e.g.
                    %   obj.method(args).field
                    % so simply grab the first output and store as
                    % intermediate
                    result = obj.(name)(args{:});
                end
            else
                error('Method %s not publically accessible', name);
            end
        end

        function varargout = subsref(obj, s)
            % Handle all subsref type usage of instance, make direct property
            % and method calls, or defer to appropriate abstract methods above
            % for indexing and dynamic property/method references
            %
            % Any indexing beyond this will be passed to the builtin
            
            % for performance, the dynamic subsref can be turned off
            if ~obj.dynamicClassActive
                [varargout{1:nargout}] = builtin('subsref', obj, s);
                return;
            end
            
            % for storing intermediate results
            result = obj;
            expandResultAsCell = false;

            switch s(1).type
                case '.'
                    % obj.name...
                    name = s(1).subs;
                    meta = metaclass(obj);
                    isProp = isprop(obj, name);
                    isMethod = ismethod(obj, name);

                    % is this a defined property?
                    if ~obj.allowDynamicPropertyOverride() && isProp
                        % is it a public property
                        [result, s] = obj.mapDeclaredPropertyAccess(name, meta, s);

                    % is it a defined method?
                    elseif ~obj.allowDynamicMethodOverride() && isMethod
                        % is it a public method?
                        [result, s, returnImmediately] = ...
                            obj.mapDeclaredMethodCall(name, meta, s, nargout);
                        if returnImmediately
                            varargout = result;
                            return;
                        end

                    else
                        % not a defined property or method
                        foundSupported = false;

                        % try it as a dynamic property
                        % is there subsequent indexing that should be passed along?
                        if length(s) > 1 
                            type = s(2).type;
                            subs = s(2).subs;
                        else
                            type = '';
                            subs = {};
                        end

                        % request the property value
                        [value, indexingApplied] = ...
                            obj.mapDynamicPropertyAccess(name, ...
                            type, subs);

                        if ~isequal(value, DynamicClass.NotSupported)
                            foundSupported = true;
                            result = value;
                            if indexingApplied && ~isempty(type)
                                % drop the next ()/{}/. indexing as it has been applied already
                                s = s(3:end);
                            else
                                % keep the next ()/{}/. indexing as it has not been applied yet
                                s = s(2:end);
                            end
                        end

                        % then try it as a dynamic method if it wasn't a dynamic property 
                        if ~foundSupported
                            dynamicMethodFn = obj.mapDynamicMethod(name);
                            if ~isequal(dynamicMethodFn, DynamicClass.NotSupported) && ...
                                ~isempty(dynamicMethodFn)

                                foundSupported = true;
                               
                                % this method is supported, figure out how to call it
                                if length(s) > 1 && strcmp(s(2).type, '()')
                                    % method call with argument list
                                    args = s(2).subs;
                                    s = s(3:end);
                                else
                                    % method call without arguments
                                    args = {};
                                    s = s(2:end);
                                end

                                % does further subsref indexing happen after this?
                                if isempty(s)
                                    % no more after this, call the method and return
                                    % this syntax calls the method and preserves
                                    % the correct nargout
                                    [varargout{1:nargout}] = dynamicMethodFn(args{:});
                                    return;
                                else
                                    % more subsref indexing happens after this, e.g.
                                    %   obj.method(args).field
                                    % capture only the first output
                                    result = dynamicMethodFn(args{:});
                                end
                            end
                        end

                        if ~foundSupported
                            % now try the static property and method access again
                            if obj.allowDynamicPropertyOverride() && isProp 
                                [result, s] = obj.mapDeclaredPropertyAccess(name, meta, s);
                                
                            elseif obj.allowDynamicMethodOverride() && isMethod
                                [result, s, returnImmediately] = ...
                                    obj.mapDeclaredMethodCall(name, meta, s);
                                if returnImmediately
                                    varargout = result;
                                    return;
                                end
                                
                            else
                                % unclear what this is
                                error('Could not find defined or dynamic property or method %s', name);
                            end
                        end
                    end

                case {'()', '{}'}
                    type = s(1).type;
                    subs = s(1).subs;

                    % gather the next .field indexing to pass along
                    if length(s) > 1
                        typeNext = s(2).type;
                        subsNext = s(2).subs;
                    else
                        typeNext = '';
                        subsNext = {}; 
                    end
                        
                    if strcmp(type, '()')
                        [result, appliedNext] = obj.parenIndex(subs, typeNext, subsNext);
                        typeStr = 'Parenthetical';
                    else
                        [result, appliedNext] = obj.cellIndex(subs, typeNext, subsNext);
                        typeStr = 'Cell array';
                        expandResultAsCell = true;
                    end

                    if ~isequal(result, DynamicClass.NotSupported)
                        if ~isempty(typeNext) && appliedNext
                            s = s(3:end);
                        else
                            s = s(2:end);
                        end
                    else
                        error('%s indexing not supported for this class', typeStr);
                    end
            end

            % if there is any remaining indexing, use subsref to handle it, which
            % may result in recursion if result is a class instance
            if ~isempty(s)
                % this syntax preserves the correct number of output arguments
                if expandResultAsCell
                    result = result{1};
                end
                [varargout{1:nargout}] = subsref(result, s);
            else
                if expandResultAsCell
                    % when using cell indexing {}, result must be expanded to varargout
                    assert(iscell(result));
                    varargout = result;
                else
                    varargout{1} = result;
                end
            end
        end

        function obj = subsasgn(obj, s, value)
            if strcmp(s(1).subs, 'dynamicClassActive')
                obj.dynamicClassActive = value;
                return;
            end
            if ~obj.dynamicClassActive
                obj = builtin('subsasgn', obj, s, value);
                return;
            end
            
            sExtra = s(2:end);

            switch s(1).type;
                case '.'
                    name = s(1).subs;
                    isProp = isprop(obj, name);
                    
                    retVal = DynamicClass.NotSupported;
                    if obj.allowDynamicPropertyOverride || ~isProp
                        % it's either not a declared property, or we're allowing 
                        % dynamicPropertyAssign to handle declared properties
                        retVal = obj.dynamicPropertyAssign(name, value, sExtra);
                    end

                    if isProp && isequal(retVal, DynamicClass.NotSupported)
                        % either we're not allowing dynamic property override
                        % or dynamicPropertyAssign doesn't want to handle it
                        % access it like a declared property
                        
                        % is this a public property?
                        meta = metaclass(obj);
                        propInfo = meta.PropertyList;
                        idx = find(strcmp(name, {propInfo.Name}));
                        assert(~isempty(idx), 'Could not find property info for %s', name);
                        if strcmp(propInfo(idx).SetAccess, 'public')
                            % yes, pass to builtin below
                            retVal = DynamicClass.NotSupported;
                        else
                            error('Property %s is not publically setable', name);
                        end
                    end

                case '()'
                    subs = s(1).subs;
                    retVal = obj.parenAssign(subs, value, sExtra);

                case '{}'
                    subs = s(1).subs;
                    retVal = obj.cellAssign(subs, value, sExtra);
            end

            if isequal(retVal, DynamicClass.NotSupported)
                obj = builtin('subsasgn', obj, s, value);
            else
                obj = retVal;
            end
        end

    end
end
