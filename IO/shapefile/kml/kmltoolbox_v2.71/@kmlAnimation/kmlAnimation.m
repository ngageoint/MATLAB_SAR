classdef kmlAnimation < handle
%KMLANIMATION(name) Create a KML animation helper, allowing you to create timed animations.
%   This class should not be instantiated directly, but instead through kml.newAnimation
%    
%   Copyright 2013 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.6 $  $Date: 2013/05/17 17:17:17 $
   
    properties
      kml 
      name
      tour
      playlist
   end
   
   methods
       function this = kmlAnimation(KML, Name)
           assert(isa(KML,'kml'),'Invalid usage of kmlAnimation');
           this.kml  = KML;
           this.name = Name;
           
           this.tour     = this.kml.xml.createElement('gx:Tour');
           this.playlist = this.kml.xml.createElement('gx:Playlist');
           this.tour.appendChild(this.kml.textNode('name',this.name));
           this.tour.appendChild(this.playlist);
           this.kml.doc.appendChild(this.tour);
       end
       
       
       function updateLocation(this,target,waitingTime,longitude,latitude,altitude)
           if numel(target)>1
              if numel(longitude)==1
                  longitude = repmat(longitude,numel(target),1);
              end
              if numel(latitude)==1
                  latitude = repmat(latitude,numel(target),1);
              end
              if numel(altitude)==1
                  altitude = repmat(altitude,numel(target),1);
              end
              for i = 1:numel(target)
                  this.updateLocation(target(i),waitingTime,longitude(i),latitude(i),altitude(i));
              end
              return
           end
           
           assert(isfield(target,'location_id'),'Invalid target for location update');
           
           kml = this.kml;
           [longitude,latitude] = kml.checkUnit(longitude,latitude);
            
           animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
           update      = kml.xml.createElement('Update');
           change      = kml.xml.createElement('Change');
           location    = kml.xml.createElement('Location');
           animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           
           location.setAttribute('targetId',target.location_id);
           location.appendChild(kml.textNode('latitude',  sprintf('%0.16g',latitude)));
           location.appendChild(kml.textNode('longitude', sprintf('%0.16g',longitude)));
           location.appendChild(kml.textNode('altitude',  sprintf('%0.16g',altitude)));
           
           change.appendChild(location);
           update.appendChild(change);
           animUpdate.appendChild(update);
           
           this.playlist.appendChild(animUpdate);
       end
       
       function updateOrientation(this,target,waitingTime,heading,tilt,roll)
           if numel(target)>1
              if numel(heading)==1
                  heading = repmat(heading,numel(target),1);
              end
              if numel(tilt)==1
                  tilt = repmat(tilt,numel(target),1);
              end
              if numel(roll)==1
                  roll = repmat(roll,numel(target),1);
              end
              for i = 1:numel(target)
                  this.updateOrientation(target(i),waitingTime,heading(i),tilt(i),roll(i));
              end
              return
           end
           
           
           assert(isfield(target,'orientation_id'),'Invalid target for orientation update')
           kml = this.kml;
           [heading,tilt,roll] = kml.checkUnit(heading,tilt,roll);
            
           animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
           update      = kml.xml.createElement('Update');
           change      = kml.xml.createElement('Change');
           orientation    = kml.xml.createElement('Orientation');
           animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           
           orientation.setAttribute('targetId',target.orientation_id);
           orientation.appendChild(kml.textNode('heading', sprintf('%0.16g',heading)));
           orientation.appendChild(kml.textNode('tilt',    sprintf('%0.16g',tilt)));
           orientation.appendChild(kml.textNode('roll',    sprintf('%0.16g',roll)));
        
           change.appendChild(orientation);
           update.appendChild(change);
           animUpdate.appendChild(update);
           
           this.playlist.appendChild(animUpdate);
       end       
       
       function updateScale(this,target,waitingTime,scaleX,scaleY,scaleZ)
           if numel(target)>1
              if numel(scaleX)==1
                  scaleX = repmat(scaleX,numel(target),1);
              end
              if numel(scaleY)==1
                  scaleY = repmat(scaleY,numel(target),1);
              end
              if numel(scaleZ)==1
                  scaleZ = repmat(scaleZ,numel(target),1);
              end
              for i = 1:numel(target)
                  this.updateScale(target(i),waitingTime,scaleX(i),scaleY(i),scaleZ(i));
              end
              return
           end
           
           assert(isfield(target,'scale_id'),'Invalid target for scale update')
           kml = this.kml;
            
           animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
           update      = kml.xml.createElement('Update');
           change      = kml.xml.createElement('Change');
           scale       = kml.xml.createElement('Scale');
           animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           
           scale.setAttribute('targetId',target.scale_id);
           scale.appendChild(kml.textNode('x', sprintf('%0.16g',scaleX)));
           scale.appendChild(kml.textNode('y', sprintf('%0.16g',scaleY)));
           scale.appendChild(kml.textNode('z', sprintf('%0.16g',scaleZ)));
        
           change.appendChild(scale);
           update.appendChild(change);
           animUpdate.appendChild(update);
           
           this.playlist.appendChild(animUpdate);
       end            
       
       function updateVisibility(this,target,visibility)
           if numel(target)>1
              if numel(visibility)==1
                  visibility = repmat(visibility,numel(target),1);
              end
              for i = 1:numel(target)
                  this.updateVisibility(target(i),visibility(i));
              end
              return
           end           
           assert(isfield(target,'id'),'Invalid target for visibility update');
           assert(isfield(target,'type'),'Invalid target for visibility update');
           
           kml = this.kml;
            
           animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
           update      = kml.xml.createElement('Update');
           change      = kml.xml.createElement('Change');
           item        = kml.xml.createElement(target.type);
%            animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           
           item.setAttribute('targetId',target.id);
           item.appendChild(kml.textNode('visibility', num2str(visibility)));
        
           change.appendChild(item);
           update.appendChild(change);
           animUpdate.appendChild(update);
           
           this.playlist.appendChild(animUpdate);
       end
       
       function updateText(this,target,text)
           if numel(target)>1
              if numel(text)==1
                  text = repmat({text},numel(target),1);
              end
              for i = 1:numel(target)
                  this.updateText(target(i),text(i));
              end
              return
           end           
           assert(isfield(target,'id'),'Invalid target for visibility update');
           assert(isfield(target,'type'),'Invalid target for visibility update');
           
           kml = this.kml;
            
           animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
           update      = kml.xml.createElement('Update');
           change      = kml.xml.createElement('Change');
           item        = kml.xml.createElement(target.type);
%            animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           
           item.setAttribute('targetId',target.id);
           item.appendChild(kml.textNode('name', text));
        
           change.appendChild(item);
           update.appendChild(change);
           animUpdate.appendChild(update);
           
           this.playlist.appendChild(animUpdate);
       end       
       
       function updateCoordinates(this,target,waitingTime,longitude,latitude,altitude)
           kml = this.kml;
                   
           [longitude,latitude] = kml.checkUnit(longitude,latitude);
           if numel(target)>1
              if numel(longitude)==1
                  longitude = repmat(longitude,numel(target),1);
                  latitude = repmat(latitude,numel(target),1);
                  altitude = repmat(altitude,numel(target),1);
              end
           end
           assert(isfield(target,'id'),'Invalid target for coordinates update');
           assert(isfield(target,'coordinates_type'),'Invalid target for coordinates update');
           assert(isfield(target,'coordinates_id'),'Invalid target for coordinates update');

   
           
           for i = 1:numel(target)
               if numel(target)==1 && numel(latitude)>1
                   coordinates  = sprintf('%0.16g,%0.16g,%0.16g ',[longitude(:) latitude(:) altitude(:)].');
               else
                   coordinates  = sprintf('%0.16g,%0.16g,%0.16g ',[longitude(i) latitude(i) altitude(i)].');   
               end

               animUpdate  = kml.xml.createElement('gx:AnimatedUpdate');
               update      = kml.xml.createElement('Update');
               change      = kml.xml.createElement('Change');
               item        = kml.xml.createElement(target(i).coordinates_type);
               animUpdate.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));

               item.setAttribute('targetId',target(i).coordinates_id);
               item.appendChild(kml.textNode('coordinates', coordinates));

               change.appendChild(item);
               update.appendChild(change);
               animUpdate.appendChild(update);

               this.playlist.appendChild(animUpdate);
           end
       end              
       
       function flyToCamera(this,waitingTime,longitude,latitude,altitude,heading,tilt,roll,altitudeMode)
           if nargin<9
               altitudeMode = 'relativeToGround';
           else
               assert(ismember(altitudeMode,{'clampToGround','relativeToGround','absolute'}));
           end

           kml = this.kml;
           
           [longitude,latitude,heading,tilt,roll] = kml.checkUnit(longitude,latitude,heading,tilt,roll);

           flyTo = kml.xml.createElement('gx:FlyTo');
           flyTo.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           flyTo.appendChild(kml.textNode('gx:flyToMode','smooth'));
           camera = kml.xml.createElement('Camera');
           camera.appendChild(kml.textNode('longitude',sprintf('%0.16g',longitude)));
           camera.appendChild(kml.textNode('latitude',sprintf('%0.16g',latitude)));
           camera.appendChild(kml.textNode('altitude',sprintf('%0.16g',altitude)));
           h = heading-180;
           
           if heading < -360
               heading = 720+heading;
           elseif h>360
               heading = heading-720;
           end
           
           camera.appendChild(kml.textNode('heading',sprintf('%0.16g',heading)));
           camera.appendChild(kml.textNode('tilt',sprintf('%0.16g',tilt)));
           camera.appendChild(kml.textNode('roll',sprintf('%0.16g',roll)));
           camera.appendChild(kml.textNode('altitudeMode',altitudeMode));
           
           flyTo.appendChild(camera);
           this.playlist.appendChild(flyTo);
       end
       
       function flyToLookAt(this,waitingTime, longitude,latitude,altitude,altitudeMode, heading, tilt)
          if nargin<6
               altitudeMode = 'relativeToGround';
          else
              assert(ismember(altitudeMode,{'clampToGround','relativeToGround','absolute'}));
          end
           
          if nargin < 7
              heading = 0;
          end
          if nargin < 8
              tilt = 0;
          end
          
           kml = this.kml;
           
           [longitude,latitude,heading,tilt] = kml.checkUnit(longitude,latitude,heading,tilt);
           
           flyTo = kml.xml.createElement('gx:FlyTo');
           flyTo.appendChild(kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           flyTo.appendChild(kml.textNode('gx:flyToMode','smooth'));
           lookAt = kml.xml.createElement('LookAt');
           lookAt.appendChild(kml.textNode('longitude',sprintf('%0.16g',longitude)));
           lookAt.appendChild(kml.textNode('latitude',sprintf('%0.16g',latitude)));
           lookAt.appendChild(kml.textNode('altitude',sprintf('%0.16g',altitude)));
           lookAt.appendChild(kml.textNode('altitudeMode',altitudeMode));
           lookAt.appendChild(kml.textNode('heading',sprintf('%0.16g',heading)));
           lookAt.appendChild(kml.textNode('tilt',sprintf('%0.16g',tilt)));
           
           flyTo.appendChild(lookAt);
           this.playlist.appendChild(flyTo);
       end
       
       function wait(this,waitingTime)
           wait  = this.kml.xml.createElement('gx:Wait');
           wait.appendChild(this.kml.textNode('gx:duration',sprintf('%0.16g',waitingTime)));
           this.playlist.appendChild(wait);
       end
       
       function pause(this)
           pause  = this.kml.xml.createElement('gx:TourControl');
           pause.appendChild(this.kml.textNode('gx:playMode','pause'));
           this.playlist.appendChild(pause);           
       end
   end
   
end