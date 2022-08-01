<?xml version="1.0" encoding="utf-8"?>

<!--
Adapted from MathWorks provided mxdom2simplehtml.xsl
This is an XSL stylesheet which converts mscript XML files into HTML.
Use the XSLT command to perform the conversion.

Copyright 1984-2019 The MathWorks, Inc.
Copyright 2022 Cristian Le
-->

<!DOCTYPE xsl:stylesheet [ <!ENTITY nbsp "&#160;"> <!ENTITY reg "&#174;"> ]>
<xsl:stylesheet
        version="1.0"
        xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
        xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd"
        exclude-result-prefixes="mwsh">
    <xsl:output method="html"
                indent="no"
                doctype-public="-//W3C//DTD HTML 4.01 Transitional//EN"/>
    <xsl:strip-space elements="mwsh:code"/>

    <xsl:variable name="title">
        <xsl:variable name="dTitle" select="//steptitle[@style='document']"/>
        <xsl:choose>
            <xsl:when test="$dTitle">
                <xsl:value-of select="$dTitle"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="mscript/m-file"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>


    <xsl:template match="mscript">
        <html>

            <!-- head -->
            <head>
                <xsl:comment>
                    This HTML was auto-generated from MATLAB code.
                    To make changes, update the MATLAB code and republish this document.
                </xsl:comment>

                <title>
                    <xsl:value-of select="$title"/>
                </title>

                <meta name="generator">
                    <xsl:attribute name="content">MATLAB
                        <xsl:value-of select="version"/>
                    </xsl:attribute>
                </meta>
                <link rel="schema.DC" href="http://purl.org/dc/elements/1.1/"/>
                <meta name="DC.date">
                    <xsl:attribute name="content">
                        <xsl:value-of select="date"/>
                    </xsl:attribute>
                </meta>
                <meta name="DC.source">
                    <xsl:attribute name="content"><xsl:value-of select="m-file"/>.m
                    </xsl:attribute>
                </meta>

                <xsl:call-template name="include-script"/>

                <meta name="chunktype" content="refpage"/>
                <meta http-equiv="Content-Script-Type" content="text/javascript"/>

            </head>

            <body id="responsive_offcanvas" scroll="no" style="overflow: hidden">

                <xsl:call-template name="header"/>

                <div class="content_container" id="content_container" itemprop="content">
                    <section id="doc_center_content">

                        <!-- Determine if the there should be an introduction section. -->
                        <xsl:variable name="hasIntro" select="count(cell[@style = 'overview'])"/>

                        <!-- If there is an introduction, display it. -->
                        <xsl:if test="$hasIntro">
                            <!-- Title and Header -->
                            <h1 itemprop="title">
                                <xsl:attribute name="id">
                                    <xsl:value-of select="translate(cell[1]/steptitle,' ','-')"/>
                                </xsl:attribute>
                                <xsl:apply-templates select="cell[1]/steptitle"/>
                            </h1>
                            <xsl:comment>introduction</xsl:comment>
                            <div class="doc_topic_desc">
                                <div class="purpose_container">
                                    <div class="switch">
                                        <a href="javascript:void(0);" id="expandAllPage"
                                           data-allexpanded="true1">collapse all in page
                                        </a>
                                    </div>
                                    <xsl:apply-templates select="cell[1]/text"/>
                                    <!-- There can be text output if there was a parse error. -->
                                    <xsl:apply-templates select="cell[1]/mcodeoutput"/>
                                </div>
                            </div>
                            <xsl:comment>/introduction</xsl:comment>
                        </xsl:if>

                        <xsl:variable name="body-cells" select="cell[not(@style = 'overview')]"/>

                        <!-- Loop over each cell -->
                        <xsl:for-each select="$body-cells">
                            <div class="ref_sect">
                                <!-- Title of cell -->
                                <h2>
                                    <xsl:attribute name="id">
                                        <xsl:value-of select="translate(steptitle,' ','-')"/>
                                    </xsl:attribute>
                                    <xsl:apply-templates select="steptitle"/>
                                </h2>
                                <!-- Contents of each cell -->
                                <xsl:apply-templates select="text"/>
                                <xsl:apply-templates select="mcode-xmlized"/>
                                <xsl:apply-templates select="mcodeoutput|img"/>
                            </div>
                        </xsl:for-each>


                    </section>
                </div>

                <xsl:call-template name="footer"/>
            </body>
        </html>
    </xsl:template>

    <xsl:template name="include-script">
        <script type="text/javascript">
            function addNext(ind, list, path) {
                if (ind >= list.length)
                    return;
                let type = list[ind].type;
                if (type === 'script')
                    addScript(ind, list, path);
                else if (type === 'link')
                    addCSSLink(ind, list, path);
                else {
                    alert(`Unknown dynamic script being added: type=${type}. Skipping.`);
                    addNext(ind + 1, list, path);
                }
            }

            function addScript(ind, list, path, currPath = path, currInd = 0) {
                let location = list[ind].location;
                let head = document.head;
                let script = document.createElement('script');
                script.src = currPath + location;
                script.onload = () => {
                    addNext(ind + 1, list, currPath);
                }
                script.onerror = (error) => {
                    if (currInd > 12) {
                        alert(`Looped too deep script=${currPath + location} ind=${currInd}`);
                        addNext(ind + 1, list, path);
                    } else
                        addScript(ind, list, path, currPath + "../", currInd + 1);
                }
                head.appendChild(script);
            }

            function addCSSLink(ind, list, path, currPath = path, currInd = 0) {
                let location = list[ind].location;
                let media = list[ind].media;
                let head = document.head;
                let link = document.createElement('link');
                link.href = currPath + location;
                link.type = "text/css";
                link.rel = "stylesheet";
                if (media !== null)
                    link.media = media;
                link.onload = () => {
                    addNext(ind + 1, list, currPath);
                }
                link.onerror = (error) => {
                    if (currInd > 10) {
                        alert(`Looped too deep link=${currPath + location} ind=${currInd}`);
                        addNext(ind + 1, list, path);
                    } else
                        addCSSLink(ind, list, path, currPath + "../", currInd + 1);
                }
                head.appendChild(link);
            }

            list = [{
                "location": "includes/product/scripts/jquery/jquery-latest.js",
                "type": "script"
            }, {
                "location": "includes/product/css/bootstrap.min.css",
                "type": "link",
                "media": null
            }, {
                "location": "includes/product/css/site6.css",
                "type": "link",
                "media": null
            }, {
                "location": "includes/product/css/site6_lg.css",
                "type": "link",
                "media": "screen and (min-width: 1200px)"
            }, {
                "location": "includes/product/css/site6_md.css",
                "type": "link",
                "media": "screen and (min-width: 992px) and (max-width: 1199px)"
            }, {
                "location": "includes/product/css/site6_sm+xs.css",
                "type": "link",
                "media": "screen and (max-width: 991px)"
            }, {
                "location": "includes/product/css/site6_sm.css",
                "type": "link",
                "media": "screen and (min-width: 768px) and (max-width: 991px)"
            }, {
                "location": "includes/product/css/site6_xs.css",
                "type": "link",
                "media": "screen and (max-width: 767px)"
            }, {
                "location": "includes/product/css/site6_offcanvas_v2.css",
                "type": "link",
                "media": null
            }, {
                "location": "includes/shared/scripts/l10n.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/docscripts.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/f1help.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/docscripts.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/mw.imageanimation.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/jquery.highlight.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/underscore-min.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/use_platform_screenshots.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/suggest.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/overload.js",
                "type": "script"
            }, {
                "location": "includes/shared/scripts/helpservices.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/productfilter.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/jquery/jquery.mobile.custom.min.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/bootstrap.min.js",
                "type": "script"
            }, {
                "location": "includes/product/scripts/global.js",
                "type": "script"
            }, {
                "location": "includes/product/css/doc_center_base.css",
                "type": "link",
                "media": null
            }, {
                "location": "includes/product/css/doc_center_installed.css",
                "type": "link",
                "media": null
            }, {
                "location": "includes/product/css/doc_center_print.css",
                "type": "link",
                "media": "print"
            }, {
                "location": "includes/shared/equationrenderer/release/MathRenderer.js",
                "type": "script"
            }];
            addNext(0, list, "");
        </script>
    </xsl:template>

    <!-- Header -->
    <xsl:template name="header">
        <div id="doc_header_spacer" class="header"></div>
    </xsl:template>

    <!-- Footer -->
    <xsl:template name="footer">
        <p class="footer">
            <xsl:value-of select="copyright"/>
            <br/>
            <a href="https://www.mathworks.com/products/matlab/">
                Published with MATLAB&reg; R<xsl:value-of select="release"/>
            </a>
            <br/>
        </p>
    </xsl:template>


    <!-- HTML Tags in text sections -->
    <xsl:template match="p">
        <p>
            <xsl:apply-templates/>
        </p>
    </xsl:template>
    <xsl:template match="ul">
        <div class="itemizedlist">
            <ul>
                <xsl:apply-templates/>
            </ul>
        </div>
    </xsl:template>
    <xsl:template match="ol">
        <div>
            <ol>
                <xsl:apply-templates/>
            </ol>
        </div>
    </xsl:template>
    <xsl:template match="li">
        <li>
            <xsl:apply-templates/>
        </li>
    </xsl:template>
    <xsl:template match="pre">
        <xsl:choose>
            <xsl:when test="@class='error'">
                <pre class="error">
                    <xsl:apply-templates/>
                </pre>
            </xsl:when>
            <xsl:otherwise>
                <pre>
                    <xsl:apply-templates/>
                </pre>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
    <xsl:template match="b">
        <b>
            <xsl:apply-templates/>
        </b>
    </xsl:template>
    <xsl:template match="i">
        <i>
            <xsl:apply-templates/>
        </i>
    </xsl:template>
    <xsl:template match="tt">
        <tt>
            <xsl:apply-templates/>
        </tt>
    </xsl:template>
    <xsl:template match="a">
        <a>
            <xsl:attribute name="href">
                <xsl:value-of select="@href"/>
            </xsl:attribute>
            <xsl:apply-templates/>
        </a>
    </xsl:template>
    <xsl:template match="html">
        <xsl:value-of select="@text" disable-output-escaping="yes"/>
    </xsl:template>
    <xsl:template match="latex"/>

    <!-- Detecting M-Code in Comments-->
    <xsl:template match="text/mcode-xmlized">
        <pre class="language-matlab">
            <xsl:apply-templates/>
            <xsl:text><!-- g162495 -->
            </xsl:text>
        </pre>
    </xsl:template>

    <!-- Code input and output -->

    <xsl:template match="mcode-xmlized">
        <pre class="codeinput">
            <xsl:apply-templates/>
            <xsl:text><!-- g162495 -->
            </xsl:text>
        </pre>
    </xsl:template>

    <xsl:template match="mcodeoutput">
        <xsl:choose>
            <xsl:when test="concat(substring(.,0,7),substring(.,string-length(.)-7,7))='&lt;html&gt;&lt;/html&gt;'">
                <xsl:value-of select="substring(.,7,string-length(.)-14)" disable-output-escaping="yes"/>
            </xsl:when>
            <xsl:otherwise>
                <pre>
                    <xsl:attribute name="class">
                        <xsl:value-of select="@class"/>
                    </xsl:attribute>
                    <xsl:apply-templates/>
                </pre>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>


    <!-- Figure and model snapshots and equations -->
    <xsl:template match="img[@class='equation']">
        <img>
            <xsl:attribute name="src">
                <xsl:value-of select="@src"/>
            </xsl:attribute>
            <xsl:attribute name="alt">
                <xsl:value-of select="@alt"/>
            </xsl:attribute>
            <xsl:if test="@scale">
                <xsl:attribute name="style">
                    <xsl:if test="@width">
                        <xsl:text>width:</xsl:text>
                        <xsl:value-of select="@width"/>
                        <xsl:text>;</xsl:text>
                    </xsl:if>
                    <xsl:if test="@height">
                        <xsl:text>height:</xsl:text>
                        <xsl:value-of select="@height"/>
                        <xsl:text>;</xsl:text>
                    </xsl:if>
                </xsl:attribute>
            </xsl:if>
        </img>
    </xsl:template>

    <xsl:template match="img">
        <img vspace="5" hspace="5">
            <xsl:attribute name="src">
                <xsl:value-of select="@src"/>
            </xsl:attribute>
            <xsl:if test="@width or @height">
                <xsl:attribute name="style">
                    <xsl:if test="@width">
                        <xsl:text>width:</xsl:text>
                        <xsl:value-of select="@width"/>
                        <xsl:text>;</xsl:text>
                    </xsl:if>
                    <xsl:if test="@height">
                        <xsl:text>height:</xsl:text>
                        <xsl:value-of select="@height"/>
                        <xsl:text>;</xsl:text>
                    </xsl:if>
                </xsl:attribute>
            </xsl:if>
            <xsl:attribute name="alt">
                <xsl:value-of select="@alt"/>
            </xsl:attribute>
            <xsl:text/>
        </img>
    </xsl:template>

    <!-- Colors for syntax-highlighted input code -->

    <xsl:template match="mwsh:code">
        <xsl:apply-templates/>
    </xsl:template>
    <xsl:template match="mwsh:keywords">
        <span class="keyword">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>
    <xsl:template match="mwsh:strings">
        <span class="string">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>
    <xsl:template match="mwsh:comments">
        <span class="comment">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>
    <xsl:template match="mwsh:unterminated_strings">
        <span class="untermstring">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>
    <xsl:template match="mwsh:system_commands">
        <span class="syscmd">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>
    <xsl:template match="mwsh:type_section">
        <span class="typesection">
            <xsl:value-of select="."/>
        </span>
    </xsl:template>


    <!-- Footer information -->

    <xsl:template match="copyright">
        <xsl:value-of select="."/>
    </xsl:template>
    <xsl:template match="revision">
        <xsl:value-of select="."/>
    </xsl:template>

    <!-- Search and replace  -->
    <!-- From http://www.xml.com/lpt/a/2002/06/05/transforming.html -->

    <xsl:template name="globalReplace">
        <xsl:param name="outputString"/>
        <xsl:param name="target"/>
        <xsl:param name="replacement"/>
        <xsl:choose>
            <xsl:when test="contains($outputString,$target)">
                <xsl:value-of select=
                                      "concat(substring-before($outputString,$target),$replacement)"/>
                <xsl:call-template name="globalReplace">
                    <xsl:with-param name="outputString"
                                    select="substring-after($outputString,$target)"/>
                    <xsl:with-param name="target" select="$target"/>
                    <xsl:with-param name="replacement"
                                    select="$replacement"/>
                </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="$outputString"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>
</xsl:stylesheet>