/* Copyright (c) 2007 Paul Bakaus (paul.bakaus@googlemail.com) and Brandon Aaron (brandon.aaron@gmail.com || http://brandonaaron.net)
 * Dual licensed under the MIT (http://www.opensource.org/licenses/mit-license.php)
 * and GPL (http://www.opensource.org/licenses/gpl-license.php) licenses.
 *
 * $LastChangedDate: 2008-02-28 05:49:55 -0500 (Thu, 28 Feb 2008) $
 * $Rev: 4841 $
 *
 * Version: @VERSION
 *
 * Requires: jQuery 1.2+
 */
(function(B){B.dimensions={version:"@VERSION"};B.each(["Height","Width"],function(D,C){B.fn["inner"+C]=function(){if(!this[0]){return }var F=C=="Height"?"Top":"Left",E=C=="Height"?"Bottom":"Right";return this.css("display")!="none"?this[0]["client"+C]:A(this,C.toLowerCase())+A(this,"padding"+F)+A(this,"padding"+E)};B.fn["outer"+C]=function(F){if(!this[0]){return }var H=C=="Height"?"Top":"Left",E=C=="Height"?"Bottom":"Right";F=B.extend({margin:false},F||{});var G=this.css("display")!="none"?this[0]["offset"+C]:A(this,C.toLowerCase())+A(this,"border"+H+"Width")+A(this,"border"+E+"Width")+A(this,"padding"+H)+A(this,"padding"+E);
return G+(F.margin?(A(this,"margin"+H)+A(this,"margin"+E)):0)}});B.each(["Left","Top"],function(D,C){B.fn["scroll"+C]=function(E){if(!this[0]){return }return E!=undefined?this.each(function(){this==window||this==document?window.scrollTo(C=="Left"?E:B(window)["scrollLeft"](),C=="Top"?E:B(window)["scrollTop"]()):this["scroll"+C]=E}):this[0]==window||this[0]==document?self[(C=="Left"?"pageXOffset":"pageYOffset")]||B.boxModel&&document.documentElement["scroll"+C]||document.body["scroll"+C]:this[0]["scroll"+C]
}});B.fn.extend({position:function(){var H=0,G=0,F=this[0],I,C,E,D;if(F){E=this.offsetParent();I=this.offset();C=E.offset();I.top-=A(F,"marginTop");I.left-=A(F,"marginLeft");C.top+=A(E,"borderTopWidth");C.left+=A(E,"borderLeftWidth");D={top:I.top-C.top,left:I.left-C.left}}return D},offsetParent:function(){var C=this[0].offsetParent;while(C&&(!/^body|html$/i.test(C.tagName)&&B.css(C,"position")=="static")){C=C.offsetParent}return B(C)}});function A(C,D){return parseInt(B.curCSS(C.jquery?C[0]:C,D,true))||0
}})(jQuery);
